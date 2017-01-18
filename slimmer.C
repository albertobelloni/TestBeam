#include <TFile.h>
#include <TFitResult.h>
#include <TF1.h>
#include <TTree.h>
#include <TGraph.h>
#include <TString.h>
#include <vector>
#include <iostream>

using std::vector;

// Correction to wire-chamber positions
const float wc_xab = +0.642158;
const float wc_xbc = +0.209902;
const float wc_xca = -0.830777;
const float wc_yab = -8.84928 ;
const float wc_ybc = +18.0467 ;
const float wc_yca = -9.20219 ;

void slimmer(const char* filename = "") {

  TFile* f_in = TFile::Open(filename);

  std::cout << "Running over " << filename << "\n";

  TTree* qie = (TTree*)f_in->Get("QIE11Data/Events");
  TTree* wc  = (TTree*)f_in->Get("WCData/Events");

  wc->AddFriend(qie);

  // Get the branches I need:
  vector<double> *xa = 0,*xb = 0,*xc = 0,*ya = 0,*yb = 0,*yc = 0;
  wc->SetBranchAddress("xA",&xa);
  wc->SetBranchAddress("xB",&xb);
  wc->SetBranchAddress("xC",&xc);
  wc->SetBranchAddress("yA",&ya);
  wc->SetBranchAddress("yB",&yb);
  wc->SetBranchAddress("yC",&yc);

  double pulse[29][50],ped[29];
  double pulse_adc[29][50],ped_adc[29];
  double pulse_tdc[29][50];
  qie->SetBranchAddress("pulse",&pulse);
  qie->SetBranchAddress("ped",&ped);
  qie->SetBranchAddress("pulse_adc",&pulse_adc);
  qie->SetBranchAddress("ped_adc",&ped_adc);
  qie->SetBranchAddress("pulse_tdc",&pulse_tdc);

  // Create new tree with only variables I need
  // I do not change their type, even if it is unconvenient, to make sure
  // trees can be used interchangeably
  vector<double>
    *xa__N = 0,*xb__N = 0,*xc__N = 0,
    *ya__N = 0,*yb__N = 0,*yc__N = 0;
  double pulse__N[29][10],ped__N[29];
  double pulse_adc__N[29][10],ped_adc__N[29];
  double pulse_tdc__N[29][10];

  double slopeX__N, interceptX__N; // new variables: tracking
  double slopeY__N, interceptY__N; // new variables: tracking

  //TFile* f_out = TFile::Open(Form("%s_slim.root",filename),"RECREATE");
  TFile* f_out = TFile::Open(TString(filename).
                             ReplaceAll(".root","_slim.root").
                             Data(),"RECREATE");
  TTree* tree = new TTree("slim","slimmed tree");

  tree->Branch("xA",&xa__N);
  tree->Branch("xB",&xb__N);
  tree->Branch("xC",&xc__N);
  tree->Branch("yA",&ya__N);
  tree->Branch("yB",&yb__N);
  tree->Branch("yC",&yc__N);

  tree->Branch("pulse",&pulse__N,"pulse[29][10]/D");
  tree->Branch("ped",&ped__N,"ped[29]/D");
  tree->Branch("pulse_adc",&pulse_adc__N,"pulse_adc[29][10]/D");
  tree->Branch("ped_adc",&ped_adc__N,"ped_adc[29]/D");
  tree->Branch("pulse_tdc",&pulse_tdc__N,"pulse_tdc[29][10]/D");

  tree->Branch("slopeX",&slopeX__N,"slopeX/D");
  tree->Branch("interceptX",&interceptX__N,"interceptX/D");
  tree->Branch("slopeY",&slopeY__N,"slopeY/D");
  tree->Branch("interceptY",&interceptY__N,"interceptY/D");

  // Some numbers we need to
  // calculate the position of the hit with Jeff's fitter
  // We will correct B and C positions, assuming A is ok
  double z_ab =  410; // distance, in mm, between A and B
  double z_ac = 2415; // distance, in mm, between A and C
  double z_ex = 7300; // distance, in mm, between A and the tiles
  double xtrack[3],ytrack[3],ztrack[3] = {0,z_ab,z_ac};

  // Now we loop on the original ntuples, and fill the new one if needed
  // (i.e., if we have only one muon per WC, and good fit of X and Y straight lines
  for (unsigned int i=0;i<qie->GetEntries();++i) {

    if (i%10000==0)
      std::cout << "  Done with " << i << " events\n";

    qie->GetEntry(i);
    wc->GetEntry(i);

    // require one and only one muon
    if (xa->size()!=1||
        xb->size()!=1||
        xc->size()!=1||
        ya->size()!=1||
        yb->size()!=1||
        yc->size()!=1)
      continue;

    // We fill the track parameters, and find the best-fit straight line
    // note that we assume uncertainty to be the same on each measurement
    xtrack[0]=xa->at(0);
    xtrack[1]=xb->at(0)-wc_xab;
    xtrack[2]=xc->at(0)+wc_xca;
    ytrack[0]=ya->at(0);
    ytrack[1]=yb->at(0)-wc_yab;
    ytrack[2]=yc->at(0)+wc_yca;

    TF1 funx("funx","pol1",0,2500);
    TGraph grx(3,ztrack,xtrack);
    grx.Fit("funx","WCQ");
    double* line_x = funx.GetParameters();
    //const double *line_x = grx.Fit("pol1","WQNCS")->GetParams();

    TF1 funy("funy","pol1",0,2500);
    TGraph gry(3,ztrack,ytrack);
    gry.Fit("funy","WCQ");
    double* line_y = funy.GetParameters();
    //const double *line_y = gry.Fit("pol1","WQNCS")->GetParams();

    // if (i%10000==0)
    //   std::cout << " xtrack: "
    //          << xtrack[0] << " "
    //          << xtrack[1] << " "
    //          << xtrack[2] << " "
    //          << " --> " << line_x[0] << " " << line_x[1] << "\n"
    //          << " ytrack: "
    //          << ytrack[0] << " "
    //          << ytrack[1] << " "
    //          << ytrack[2] << " "
    //                  << " --> " << line_y[0] << " " << line_y[1] << "\n";

    // All good, let's set the variables for the new tree

    interceptX__N = line_x[0];
    interceptY__N = line_y[0];
    slopeX__N = line_x[1];
    slopeY__N = line_y[1];

    xa__N = xa;
    xb__N = xb;
    xc__N = xc;
    ya__N = ya;
    yb__N = yb;
    yc__N = yc;

    // more painful for plain arrays
    for (int j=0;j<29;++j) {
      ped__N[j] = ped[j];
      ped_adc__N[j] = ped_adc[j];
      for (int k=0;k<10;++k) {
        pulse__N[j][k] = pulse[j][k];
        pulse_adc__N[j][k] = pulse_adc[j][k];
        pulse_tdc__N[j][k] = pulse_tdc[j][k];
      }
    }

    tree->Fill();
  }

  tree->Write(0,TObject::kOverwrite);

}
