/*/////////////////////////////////////////////////////////////////////////////
DOCUMENTATION
Twiki with installation instructions and how to run the RooFit fitter:
- https://twiki.cern.ch/twiki/bin/view/CMS/UMDHGCalSiPMCalib

Quick summary:
    export SCRAM_ARCH=slc6_amd64_gcc700  
    cmsrel CMSSW_10_3_1_patch1
    cd CMSSW_10_3_1_patch1/src 
    git clone https://github.com/yimuchen/UserUtils 
    git clone https://github.com/yimuchen/SiPMCalib
    scram b -j 12

Typical usage (unbinned fit, fit function w/o afterpulsing nor dark current):

root -l
root [1] .L fitter2017.C++
root [2] fitter2017("energy_hists.root","energy_tree_EJ_260",0,1)

A binned fit with rebinning (4 bins merge into 1), using fit model w/ afterpulsing and dark current:

root -l 
root [1] .L fitter2017.C++
root [2] fitter2017("energy_hists.root","en_bins_EJ_260",2,0,4)

 TFile*         energy_hists.root
  KEY: TH1F     en_bins_EJ_260;1
  KEY: TH1F     en_bins_EJ_260_2P;1
  KEY: TH1F     en_bins_EJ_200;1
  KEY: TH1F     en_bins_SCSN_81F1;1
  KEY: TH1F     en_bins_SCSN_81F2;1
  KEY: TH1F     en_bins_SCSN_81F3;1
  KEY: TH1F     en_bins_SCSN_81F4;1
  KEY: TH1F     en_bins_SCSN_81S;1

  KEY: TTree    energy_tree_EJ_260;1
  KEY: TTree    energy_tree_EJ_260_2P;1
  KEY: TTree    energy_tree_EJ_200;1
  KEY: TTree    energy_tree_SCSN_81F1;1
  KEY: TTree    energy_tree_SCSN_81F2;1
  KEY: TTree    energy_tree_SCSN_81F3;1
  KEY: TTree    energy_tree_SCSN_81F4;1
  KEY: TTree    energy_tree_SCSN_81S;1

 */////////////////////////////////////////////////////////////////////////////

#include "RooRealVar.h"
#include "TSystem.h"
#include "TH1F.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooWorkspace.h"

#include "SiPMCalib/SiPMCalc/interface/SiPMPdf.hpp"
#include <cstdlib>

// FULL fit includes afterpulsing and dark current parameters
enum {BASELINE=0, WITHAFTERPULSING, FULL};
enum {BINNED=0, UNBINNED};

void fitter2017(const char* filename, const char* histname,
		int fitmodel=BASELINE, int fittype=UNBINNED,
		const int rebin=1) {

  // This might need to be loaded before calling the macro file.
  // If you are running a script without the compiling flag
  const std::string libfile
    = std::getenv("CMSSW_BASE") + std::string("/lib/")
      + std::getenv("SCRAM_ARCH") + std::string("/libSiPMCalibSiPMCalc.so");
  gSystem->Load( libfile.c_str() );

  // Declaring the Variables.
  RooRealVar x( "x", "x", -50, 600 );
  RooRealVar ped( "ped", "ped",            -2,  -10, 10 );
  RooRealVar gain( "gain", "gain",         40,    0, 1000 );
  RooRealVar s0( "s0", "s0",               12,  0.1, 100 );
  RooRealVar s1( "s1", "s1",             0.01,    0, 1 );
  RooRealVar mean( "mean", "mean",         25, 0.01, 50 );
  RooRealVar lambda( "lambda", "lambda",  0.4,    0, 0.50 );

  RooRealVar alpha( "alpha", "alpha", 0, 0.5 );
  RooRealVar beta ( "beta", "beta", 10, 20000 );
  RooRealVar dcfrac( "dcfrac", "dcfrac", 0, 0.4 );
  RooRealVar eps( "eps", "eps", 1e-5, 1e-1 );

  SiPMPdf *p;

  if (fitmodel==BASELINE)
    p = new SiPMPdf( "p", "p", x, ped, gain, s0, s1,mean, lambda);
  else if (fitmodel==WITHAFTERPULSING)
    p = new SiPMPdf( "p", "p", x, ped, gain, s0, s1,mean, lambda,
		     alpha, beta);
  else if (fitmodel==FULL)
    p = new SiPMPdf( "p", "p", x, ped, gain, s0, s1,mean, lambda,
		     alpha, beta, dcfrac, eps);
  else {
    std::cout << "The SiPMPdf pointer is null: what is going on?\n";
    return;
  }

  RooWorkspace *workspace = new RooWorkspace("w","fitter2017_workspace");

  if (fittype==BINNED) {

    // OPEN HISTOGRAM
    TFile *file = TFile::Open(filename);
    TH1F *hist = (TH1F*)file->Get(histname);
    if(!hist) {
      std::cout << "The histogram pointer is null:"
	" is the histogram name correct?\n";
      return;
    }
    hist->Rebin(rebin);

    if(hist->GetDimension()==0) {
      std::cout << "The dimension of the histogram is 0:"
	" are you sure you are not passing a tree pointer here?\n";
      return;
    }

    RooDataHist data("data","data", RooArgSet(x), hist );

    // Running fit independent estimates on ped, gain, s0, s1, mean, lambda
    p->RunEstimate( data ) ;
    p->fitTo( data) ;       // Running the fit.

    RooPlot *frame = x.frame();
    data.plotOn(frame);
    p->plotOn(frame);
    frame->Draw();

    workspace->import(data);
    workspace->import(*p);

    int npars = 6 + 2*fitmodel;
    std::cout << std::endl
	      << "CHI2/NDF:  " << frame->chiSquare(npars) << std::endl
	      << "CHI2 Prob: " << TMath::Prob(frame->chiSquare()*
					      hist->GetNbinsX(),
					      hist->GetNbinsX()-npars)
	      << std::endl;
  }

  if (fittype==UNBINNED) {

    // OPEN TREE
    TFile *file = TFile::Open(filename);
    TTree *tree = (TTree*)file->Get(histname);
    if(!tree) {
      std::cout << "The tree pointer is null: is the tree name correct?\n";
      return;
    }
    RooDataSet data("data","data", tree, RooArgSet(x), "x>-50 && x<600" );

    // Running fit independent estimates on ped, gain, s0, s1, mean, lambda
    p->RunEstimate( data ) ;
    p->fitTo( data) ;       // Running the fit.

    RooPlot *frame = x.frame();
    data.plotOn(frame);
    p->plotOn(frame);
    frame->Draw();

    workspace->import(data);
    workspace->import(*p);
  }

  workspace->writeToFile(Form("roofit_%s_%d_%d_%d.root",
			      histname,fitmodel,fittype,rebin));

}
