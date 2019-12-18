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

R__LOAD_LIBRARY(libSiPMCalibSiPMCalc)
R__LOAD_LIBRARY(libSiPMCalibInvSqCalc)

/*/////////////////////////////////////////////////////////////////////////////

Typical usage:

[] .L make_fitter_plots.C++
[] make_fitter_plots("roofit_en_bins_SCSN_81F1_0_0_5.root")

The roofit file shall contain a workspace with data and fit parameters
(fit performed with SiPMPdf object)

PROBLEMS:
- vane attempt at changing X axis range of RooPlot frame object: it does
modify the axis titles and offset, just seems not to like the range change

 */////////////////////////////////////////////////////////////////////////////


void make_fitter_plots(const char* filename) {

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

  frame->GetXaxis()->SetTitle("Energy [fC]");
  frame->GetXaxis()->SetTitleOffset(1.0);
  if(!strchr(filename,'F')) // x range for finger tiles
    frame->GetXaxis()->SetRangeUser(-50,350);
  else // x range for sigma tiles
    frame->GetXaxis()->SetRangeUser(-50,550);
  frame->GetYaxis()->SetTitle("Entries");

  TCanvas* canvas = new TCanvas("canvas","",800,500);
  frame->Draw();
  canvas->Print(TString(filename).
		ReplaceAll("roofit","canvas").ReplaceAll(".root",".png"));
  canvas->Print(TString(filename).
		ReplaceAll("roofit","canvas").ReplaceAll(".root",".pdf"));
  // printing to .C does not work well; RooFit automatically
  // adds "[x]" to object names... (the last ReplaceAll is left for consistency
  canvas->Print(TString(filename).
		ReplaceAll("roofit","canvas").ReplaceAll(".root",".root"));

  std::cout << "\n\n"
	    << " Tile/fit: " << filename << "\n\n"
	    << " MEAN NUMBER P.E.: " << mean->getValV() 
	    << " +/- " << mean->getError()
	    << "\n\n";


  return;

}
