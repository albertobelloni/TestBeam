/*//////////////////////////////////////////////////////////////////////////////
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

 *///////////////////////////////////////////////////////////////////////////////

#include "RooRealVar.h"
#include "TSystem.h"
#include "TH1F.h"
#include "TFile.h"
#include "TMath.h"
#include "RooDataHist.h"
#include "RooPlot.h"

#include "SiPMCalib/SiPMCalc/interface/SiPMPdf.hpp"
#include <cstdlib>

void fitter2017(const char* filename, const char* histname, const int rebin=1) {

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
  //RooRealVar alpha( "alpha", "alpha", 0, 0.5 );
  //RooRealVar beta( "beta", "beta", 10, 20000 );
  //RooRealVar dcfrac( "dcfrac", "dcfrac", 0, 0.4 );
  //RooRealVar eps( "eps", "eps", 1e-5, 1e-1 );

  SiPMPdf p( "p", "p", x, ped, gain, s0, s1,
             mean, lambda//,
             //alpha, beta,
             //dcfrac, eps  
  );

  TFile *file = TFile::Open(filename);
  TH1F *hist = (TH1F*)file->Get(histname);
  if(!hist) {
    std::cout << "The histogram pointer is null: is the histogram name correct?\n";
    return;
  }

  hist->Rebin(rebin);

  RooDataHist data("data","data", RooArgSet(x), hist );

  p.RunEstimate( data ) ; // Running fit independent estimates on ped, gain, s0, s1, mean, lambda  
  p.fitTo( data) ;       // Running the fit.

  RooPlot *frame = x.frame();
  data.plotOn(frame);
  p.plotOn(frame);
  frame->Draw();
  std::cout << std::endl 
	    << "CHI2/NDF:  " << frame->chiSquare(6) << std::endl
	    << "CHI2 Prob: " << TMath::Prob(frame->chiSquare()*hist->GetNbinsX(),hist->GetNbinsX()-6) << std::endl;

}
