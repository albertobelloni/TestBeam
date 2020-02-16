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
  frame->GetXaxis()->SetTitleOffset(1.0);
  if(!strchr(filename,'F')) // x range for finger tiles
    frame->GetXaxis()->SetRangeUser(-50,350);
  else // x range for sigma tiles
    frame->GetXaxis()->SetRangeUser(-50,550);
  frame->GetYaxis()->SetTitle("Entries");

  TCanvas* canvas = new TCanvas("canvas","",800,500);
  canvas->SetLogy(logy);
  frame->Draw();
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

  std::cout << "\n\n"
	    << " Tile/fit: " << filename << "\n"
	    << "\n MEAN NUMBER P.E.: " << mean->getValV()
	    << " +/- " << mean->getError()
	    << "\n CROSS-TALK PROBABILITY: " << lambda->getValV()
	    << " +/- " << lambda->getError()
	    << "\n\n";

  // Save formatted outputs to a file
  FILE *output;
  output = fopen("yield_results.txt","a");

  // Formatted output for presentation
  int row=idx(filename);
  fprintf(output,"PRESENTATION table.cell(%d,1).text = \"%.3f+-%.3f\"\n",row,
	 mean->getValV(),mean->getError());
  fprintf(output,"PRESENTATION table.cell(%d,2).text = \"%.3f+-%.3f\"\n",row,
	 lambda->getValV(),lambda->getError());
  fprintf(output,"PRESENTATION table.cell(%d,3).text = \"%.3f+-%.3f\"\n",row,
	  mean->getValV()/(1-lambda->getValV()),
	  mean->getValV()/(1-lambda->getValV())*
	  sqrt((mean->getError()/mean->getValV())*
	       (mean->getError()/mean->getValV())+
	       (lambda->getError()/(1-lambda->getValV())*
		lambda->getError()/(1-lambda->getValV()))));

  // Formatted output for DN-18-007
  fprintf(output,
	  "DN-18-007_NOTE %d &$%.3f\\pm%.3f$&$%.3f\\pm%.3f$&$%.3f\\pm%.3f$\n",
	  row,
	  mean->getValV(),mean->getError(),
	  lambda->getValV(),lambda->getError(),
	  mean->getValV()/(1-lambda->getValV()),
	  mean->getValV()/(1-lambda->getValV())*
	  sqrt((mean->getError()/mean->getValV())*
	       (mean->getError()/mean->getValV())+
	       (lambda->getError()/(1-lambda->getValV())*
		lambda->getError()/(1-lambda->getValV()))));

  fclose(output);

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

  for (int i=0;i<8;make_fitter_plots(filenames[i++]));

  return;

}

int idx(const char* filename) {

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
