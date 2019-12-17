#include "RooRealVar.h"
#include "TSystem.h"
#include "TH1F.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooPlot.h"

#include "SiPMCalib/SiPMCalc/interface/SiPMPdf.hpp"
#include <cstdlib>
#include "fitter2017.C"

R__LOAD_LIBRARY(libSiPMCalibSiPMCalc)
R__LOAD_LIBRARY(libSiPMCalibInvSqCalc)


void run_fitter(const char* tile = "energy_tree_SCSN_81S") {
  TDatime start;
  start.Print();
  gSystem->Load("fitter2017_C.so");
  fitter2017("energy_hists.root",tile);
  TDatime stop;
  stop.Print();

}
