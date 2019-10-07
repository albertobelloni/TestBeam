#include <TFile.h>
#include <TFitResult.h>
#include <TF1.h>
#include <TTree.h>
#include <TGraph.h>
#include <TString.h>
#include <vector>
#include <iostream>
#include <string>

using std::vector;

// Correction to wire-chamber positions
const double wc_xca = -1.21427;
const double wc_yca = -10.0635;

// not used, but why not saving these
const double wc_xca_res = 1.97362;
const double wc_yca_res = 1.06969;

// a couple of magic numbers
const int NCH = 16;
const int NSL = 10; // I think this is the number of time slices: in 2017 we had 50, probably because TDC had 50 (but ADC only 10)


void slimmer_2017(const char* filename = "") {

  TFile* f_in = TFile::Open(filename);

  std::cout << "Running over " << filename << "\n";

  TTree* qie = (TTree*)f_in->Get("QIE11Data/Events");
  TTree* wc  = (TTree*)f_in->Get("WCData/Events");

  wc->AddFriend(qie);

  // Get the branches I need:
  vector<double> *xa = 0,*xc = 0,*ya = 0,*yb = 0,*yc = 0;
  wc->SetBranchAddress("xA",&xa);
  wc->SetBranchAddress("xC",&xc);
  wc->SetBranchAddress("yA",&ya);
  wc->SetBranchAddress("yC",&yc);

  float pulse[NCH][NSL],ped[NCH];
  unsigned char pulse_adc[NCH][NSL];
  float ped_adc[NCH];

  qie->SetBranchAddress("pulse",&pulse);
  qie->SetBranchAddress("ped",&ped);
  qie->SetBranchAddress("pulse_adc",&pulse_adc);
  qie->SetBranchAddress("ped_adc",&ped_adc);

  // Create new tree with only variables I need
  // I do not change their type, even if it is unconvenient, to make sure
  // trees can be used interchangeably
  vector<double>
    *xa__N = 0,           *xc__N = 0,
    *ya__N = 0,           *yc__N = 0;
  float pulse__N[NCH][NSL],ped__N[NCH];
  float pulse_adc__N[NCH][NSL],ped_adc__N[NCH];

  // JAMES N ADDED TO ADD NAME TO TREE
  string run = filename;
  run = run.substr(run.find("run00") + 5, 4);
  int run__N = TString(run).Atoi();
  std::cout << "Run number: " << run__N << std::endl;

  double slopeX__N, interceptX__N; // new variables: tracking
  double slopeY__N, interceptY__N; // new variables: tracking

  TFile* f_out = TFile::Open(TString(filename).
                             ReplaceAll(".root","_slim.root").
                             Data(),"RECREATE");
  TTree* tree = new TTree("slim","slimmed tree");

  // Add a branch to store the run number
  tree->Branch("run", &run__N,"run/I");

  tree->Branch("xA",&xa__N);
  tree->Branch("xC",&xc__N);
  tree->Branch("yA",&ya__N);
  tree->Branch("yC",&yc__N);

  tree->Branch("pulse",&pulse__N,Form("pulse[%i][%i]/F",NCH,NSL));
  tree->Branch("ped",&ped__N,Form("ped[%i]/F",NCH));
  tree->Branch("pulse_adc",&pulse_adc__N,Form("pulse_adc[%i][%i]/F",NCH,NSL));
  tree->Branch("ped_adc",&ped_adc__N,Form("ped_adc[%i]/F",NCH));

  tree->Branch("slopeX",&slopeX__N,"slopeX/D");
  tree->Branch("interceptX",&interceptX__N,"interceptX/D");
  tree->Branch("slopeY",&slopeY__N,"slopeY/D");
  tree->Branch("interceptY",&interceptY__N,"interceptY/D");

  // Some numbers we need to calculate the position of the hit with Jeff's fitter
  // We will correct B and C positions, assuming A is ok
  // There is no B wire chamber in 2017 data, fit is trivial: a line through two points...
  //double z_ab =  410; // distance, in mm, between A and B; kept for bookkeeping
  double z_ac = 2415; // distance, in mm, between A and C
  double z_ex = 7300; // distance, in mm, between A and the tiles
  double xtrack[2],ytrack[2],ztrack[2] = {0,z_ac};

  // Now we loop on the original ntuples, and fill the new one if needed
  // (i.e., if we have only one muon per WC, and good fit of X and Y straight lines
  for (unsigned int i=0;i<qie->GetEntries();++i) {

    if (i%10000==0)
      std::cout << "  Done with " << i << " events\n";

    qie->GetEntry(i);
    wc->GetEntry(i);

    // require one and only one muon
    // using only coordinate y in wire chamber B
    if (xa->size()!=1||
        xc->size()!=1||
        ya->size()!=1||
        yc->size()!=1)
      continue;

    // We fill the track parameters, and find the best-fit straight line
    // note that we assume uncertainty to be the same on each measurement
    xtrack[0]=xa->at(0);
    xtrack[1]=xc->at(0)+wc_xca;
    ytrack[0]=ya->at(0);
    ytrack[1]=yc->at(0)+wc_yca;

    TF1 funx("funx","pol1",0,2500);
    TGraph grx(2,ztrack,xtrack);
    grx.Fit("funx","WCQ");
    double* line_x = funx.GetParameters();

    TF1 funy("funy","pol1",0,2500);
    TGraph gry(2,ztrack,ytrack);
    gry.Fit("funy","WCQ");
    double* line_y = funy.GetParameters();

    interceptX__N = line_x[0];
    interceptY__N = line_y[0];
    slopeX__N = line_x[1];
    slopeY__N = line_y[1];

    xa__N = xa;
    xc__N = xc;
    ya__N = ya;
    yc__N = yc;

    // more painful for plain arrays
    for (int j=0;j<NCH;++j) {
      ped__N[j] = ped[j];
      ped_adc__N[j] = ped_adc[j];
      for (int k=0;k<NSL;++k) {
        pulse__N[j][k] = pulse[j][k];
        pulse_adc__N[j][k] = pulse_adc[j][k];
      }
    }

    tree->Fill();
  }

  tree->Write(0,TObject::kOverwrite);

}
