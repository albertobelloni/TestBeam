/*/////////////////////////////////////////////////////////////////////////////

Typical usage:

[] .L make_fitter_plots.C++
[] make_fitter_plots("roofit_en_bins_SCSN_81F1_0_0_5.root")

One can also produce all plots with:

[] .L make_fitter_plots.C++
[] make_all_fitter_plots()

The roofit file shall contain a workspace with data and fit parameters
(fit performed with SiPMPdf object)

PROBLEMS:
- vane attempt at changing X axis range of RooPlot frame object: it does
modify the axis titles and offset, just seems not to like the range change

 */////////////////////////////////////////////////////////////////////////////

#include "RooRealVar.h"
#include "TSystem.h"
#include "TH1F.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooWorkspace.h"

#include "SiPMCalib/SiPMCalc/interface/SiPMPdf.hpp"
#include <cstdlib>
#include <cstring>
#include <stdio.h>

R__LOAD_LIBRARY(libSiPMCalibSiPMCalc)
R__LOAD_LIBRARY(libSiPMCalibInvSqCalc)

// Helper function used to print the fit results
// for the DN-18-007 note and the HCAL DPG presentation
int idx(const char* treename);
const char* tile(const char* filename);

void make_fitter_plots(const char* filename, bool logy = false) {

  TFile *file = TFile::Open(filename);
  if (!file) {
    std::cout << "File " << filename << " does not seem to exist\n";
    return;
  }

  // This assumes that the workspace name is "w"
  RooWorkspace *workspace = (RooWorkspace*)file->Get("w");

  RooRealVar *x      = workspace->var("x");

  RooRealVar *ped    = workspace->var("ped");
  RooRealVar *gain   = workspace->var("gain");
  RooRealVar *s0     = workspace->var("s0");
  RooRealVar *s1     = workspace->var("s1");
  RooRealVar *mean   = workspace->var("mean");
  RooRealVar *lambda = workspace->var("lambda");
  RooRealVar *alpha  = workspace->var("alpha");
  RooRealVar *beta_  = workspace->var("beta"); // beta conflicts with std::beta
  RooRealVar *dcfrac = workspace->var("dcfrac");
  RooRealVar *eps    = workspace->var("eps");


  SiPMPdf *pdf;
  if(!alpha)
    pdf = new SiPMPdf( "p", "p", *x, *ped, *gain, *s0, *s1, *mean, *lambda);
  else if (!dcfrac)
    pdf = new SiPMPdf( "p", "p", *x, *ped, *gain, *s0, *s1, *mean, *lambda, 
		       *alpha, *beta_);
  else
    pdf = new SiPMPdf( "p", "p", *x, *ped, *gain, *s0, *s1, *mean, *lambda,
		       *alpha, *beta_, *dcfrac, *eps);
  
  RooAbsData *data = workspace->data("data");

  RooPlot *frame = x->frame();

  data->plotOn(frame);
  pdf->plotOn(frame);

  frame->GetXaxis()->SetTitle("Charge [fC]");
  frame->GetXaxis()->SetTitleOffset(1.1);

  frame->GetXaxis()->SetRangeUser(-50,400);
  // Change the X range only if finger-tile linear-scale plot
  if (!logy && strchr(filename,'F'))
      frame->GetXaxis()->SetRangeUser(-50,350);

  frame->GetYaxis()->SetTitle("Events");

  frame->GetYaxis()->SetTitleOffset(1.0);
  frame->GetYaxis()->SetMaxDigits(3);

  // Let us calculate here some numbers, which we will print in the summary text file
  // and probably write also on the histogram
  float meanFit = mean->getValV();
  float meanFitErr = mean->getError();
  float lambdaFit = lambda->getValV();
  float lambdaFitErr = lambda->getError();

  float meanFitCorr = meanFit/(1-lambdaFit);
  float meanFitCorrErr = meanFit/(1-lambdaFit)*
    sqrt((meanFitErr/meanFit)*
	 (meanFitErr/meanFit)+
	 (lambdaFitErr/(1-lambdaFit)*
	  lambdaFitErr/(1-lambdaFit)));

  TCanvas* canvas = new TCanvas("canvas","",800,600);
  canvas->SetLogy(logy);
  frame->Draw();

  // Started with idea of writing also mean and lambda on plot, but there is not much space...
  TLatex label;
  label.SetNDC();
  label.SetTextSize(0.05);
  label.SetTextAlign(31); // right-aligned; center-aligned
  //label.DrawLatex(0.91,0.875,Form("mean #mu: %.3f#pm%.3f",meanFit,meanFitErr));
  //label.DrawLatex(0.91,0.825,Form("cross-talk #chi: %.3f#pm%.3f",lambdaFit,lambdaFitErr));
  //label.DrawLatex(0.91,0.775,Form("#mu/(1-#chi): %.3f#pm%.3f",meanFitCorr,meanFitCorrErr));
  label.DrawLatex(0.91,0.875,TString(tile(filename)).ReplaceAll("_","-"));

  canvas->Print(TString(filename).
		ReplaceAll("roofit","canvas").
		ReplaceAll(".root",Form("%s.png",logy?"_log":"_lin")));
  canvas->Print(TString(filename).
		ReplaceAll("roofit","canvas").
		ReplaceAll(".root",Form("%s.pdf",logy?"_log":"_lin")));
  // printing to .C does not work well; RooFit automatically
  // adds "[x]" to object names...
  canvas->Print(TString(filename).
		ReplaceAll("roofit","canvas").
		ReplaceAll(".root",Form("%s.root",logy?"_log":"_lin")));

  // Let us print the relevant fit results to the screen
  std::cout << "\n\n"
	    << " Tile/fit: " << filename << "\n"
	    << "\n MEAN NUMBER P.E.: " << meanFit
	    << " +/- " << meanFitErr
	    << "\n CROSS-TALK PROBABILITY: " << lambdaFit
	    << " +/- " << lambdaFitErr
	    << "\n\n";

  // Save formatted outputs to a file
  FILE *output;
  output = fopen("yield_results.txt","a");

  // Formatted output for presentation
  int row=idx(filename);
  fprintf(output,"PRESENTATION table.cell(%d,1).text = \"%.3f+-%.3f\"\n",row,
	 meanFit,meanFitErr);
  fprintf(output,"PRESENTATION table.cell(%d,2).text = \"%.3f+-%.3f\"\n",row,
	 lambdaFit,lambdaFitErr);
  fprintf(output,"PRESENTATION table.cell(%d,3).text = \"%.3f+-%.3f\"\n",row,
	  meanFitCorr, meanFitCorrErr);

  // Formatted output for DN-18-007
  fprintf(output,
	  "DN-18-007_NOTE %d &$%.3f\\pm%.3f$&$%.3f\\pm%.3f$&$%.3f\\pm%.3f$\n",
	  row, meanFit,meanFitErr, lambdaFit,lambdaFitErr, meanFitCorr,meanFitCorrErr);

  fclose(output);

  delete canvas;

  return;

}

void make_all_fitter_plots() {

  const char* filenames[8] =
    {"results/roofit_energy_tree_EJ_200_0_1_1.root",
     "results/roofit_energy_tree_EJ_260_0_1_1.root",
     "results/roofit_energy_tree_EJ_260_2P_0_1_1.root",
     "results/roofit_energy_tree_SCSN_81F1_0_1_1.root",
     "results/roofit_energy_tree_SCSN_81F2_0_1_1.root",
     "results/roofit_energy_tree_SCSN_81F3_0_1_1.root",
     "results/roofit_energy_tree_SCSN_81F4_0_1_1.root",
     "results/roofit_energy_tree_SCSN_81S_0_1_1.root"};

  // First make all plots in linear scale, then in logarithmic scale
  for (int i=0;i<8;make_fitter_plots(filenames[i++],false));
  for (int i=0;i<8;make_fitter_plots(filenames[i++],true));

  return;

}

int idx(const char* filename) {

  // We count tiles from 1 to 8 here, because we use idx to select
  // the row number in the yield table in the pptx presentation

  const char* tilenames[8] =
    {"EJ_200",
     "EJ_260",
     "EJ_260_2P",
     "SCSN_81F1",
     "SCSN_81F2",
     "SCSN_81F3",
     "SCSN_81F4",
     "SCSN_81S"};

  int i = 0;
  for (i=0;i<8;++i)
    if (strstr(filename,tilenames[i])!=0) break;

  // Argh, EJ_260 matches also the EJ_260_2P file, let us fix this
  // by returning immediately the row corresponding to EJ_260_2P
  if (strstr(filename,"EJ_260_2P")!=0) return 3;

  return i+1;

}

const char* tile(const char* filename) {

  const char* tilenames[8] =
    {"EJ_200",
     "EJ_260",
     "EJ_260_2P",
     "SCSN_81F1",
     "SCSN_81F2",
     "SCSN_81F3",
     "SCSN_81F4",
     "SCSN_81S"};

  int i = 0;
  for (i=0;i<8;++i)
    if (strstr(filename,tilenames[i])!=0) break;

  // Argh, EJ_260 matches also the EJ_260_2P file, let us fix this
  // by returning immediately the row corresponding to EJ_260_2P
  // Also note the removal of the underscore in front of 2P: I am later
  // going to replace all underscores with dashes, but I still want
  // a space between 260 and 2P
  if (strstr(filename,"EJ_260_2P")!=0) return "EJ_260 2P";

  return tilenames[i];

}
