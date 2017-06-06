// Corrections to get alignment:
// # Wire chamber means and standard deviations (xA-xC, xA-xC, yA-yC, etc.)
// wc_res = {}
// wc_res["x", "BC", "mean"] = -5.83e-01
// wc_res["y", "BC", "mean"] = -1.75e+01
// wc_res["x", "AC", "mean"] = -1.24e+00
// wc_res["y", "AC", "mean"] = -8.78e+00
// wc_res["x", "BC", "rms" ] =  3.96e+00
// wc_res["y", "BC", "rms" ] =  3.88e+00
// wc_res["x", "AC", "rms" ] =  4.30e+00
// wc_res["y", "AC", "rms" ] =  5.08e+00

// Based on numbers above, I get:
//     "xA-xB+0.657",
//     "xB-xC+0.583",
//     "xC-xA-1.24",
//     "yA-yB-1.12",
//     "yB-yC-3.96",
//     "yC-yA+5.08"

// Correction to wire-chamber positions
const float wc_xab = +0.642158;
const float wc_xbc = +0.209902;
const float wc_xca = -0.830777;
const float wc_yab = -8.84928 ;
const float wc_ybc = +18.0467 ;
const float wc_yca = -9.20219 ;

// (Possibly) using sigma of delta-position to clean up more
const float wc_xab_res = 1.16192 ;
const float wc_xbc_res = 1.05491 ;
const float wc_xca_res = 1.54024 ;
const float wc_yab_res = 1.2041  ;
const float wc_ybc_res = 0.892824;
const float wc_yca_res = 1.33447 ;

// Fiducial region: for each tile, 4 points are needed (8 numbers)
// TB: top, bottom; LR: left, right
// x,y BL; x,y BR; x,y TL; x,y TR
const float fiducialX[5][4] = {
  {-42,35,-37,38},
  {-37,40,-37,40},
  {-50,29,-49,30},
  {-50,28,-50,28},
  {-52,27,-51,28}
};
const float fiducialY[5][4] = {
  {-48,-50,25,24},
  {-45,-45,33,33},
  {-45,-47,34,32},
  {-49,-49,32,32},
  {-41,-42,36,35}
};

// Need to debug this!
bool isFiducial(int tile, float x, float y);

// distance from wire-chamber A to plastic tiles, in [mm]
const double z_ex = 7300;

void doAlignmentPlots(const char*
                      dir="/data/users/abelloni/CERN_TB_Aug15_slim_ntuples",
                      bool debug=false) {

  TChain* chain = new TChain("slim");
  chain->Add(Form("%s/*_slim.root",dir));

  // Plots without correction
  TH1F* hist[12];
  const string var[12] = {
    "xA-xB",
    "xB-xC",
    "xC-xA",
    "yA-yB",
    "yB-yC",
    "yC-yA",
    Form("xA-xB%s%f",wc_xab>0?"+":"",wc_xab),
    Form("xB-xC%s%f",wc_xbc>0?"+":"",wc_xbc),
    Form("xC-xA%s%f",wc_xca>0?"+":"",wc_xca),
    Form("yA-yB%s%f",wc_yab>0?"+":"",wc_yab),
    Form("yB-yC%s%f",wc_ybc>0?"+":"",wc_ybc),
    Form("yC-yA%s%f",wc_yca>0?"+":"",wc_yca)
  };
  const int color[12] = {
    kRed,
    kBlue,
    kGreen,
    kRed,
    kBlue,
    kGreen,
    kRed,
    kBlue,
    kGreen,
    kRed,
    kBlue,
    kGreen
  };

  double mean[12], sigma[12];
  for (int i=0;i<12;++i) {
    hist[i] = new TH1F(var[i].c_str(),"",240,-30,30);
    chain->Project(var[i].c_str(),var[i].c_str(),
                   Form("abs(%s)<30",var[i].c_str()));
    mean[i] = hist[i]->Fit("gaus","SQN")->GetParams()[1];
    // really not best way: why re-doing fit? I just need to
    // check if resolution of wire chambers are similar or not
    // (if different, will need to include uncertainties in tracking fit
    sigma[i] = hist[i]->Fit("gaus","SQN")->GetParams()[2];

    hist[i]->SetLineWidth(2);
    hist[i]->SetLineColor(color[i]);
    if (i<3||(i>5&&i<9))
      hist[i]->GetXaxis()->SetTitle("#deltax [mm]");
    else
      hist[i]->GetXaxis()->SetTitle("#deltay [mm]");
    hist[i]->GetXaxis()->SetTitleOffset(1.2);
    hist[i]->GetYaxis()->SetTitle("Events");
    hist[i]->GetYaxis()->SetTitleOffset(1.6);
  }

  for (int i=0;i<12;++i)
    cout << var[i].c_str() << " mean  -> " << mean[i] << endl;
  for (int i=0;i<12;++i)
    cout << var[i].c_str() << " sigma -> " << sigma[i] << endl;

  // Make canvases, we want 4 plots
  TCanvas* canv[4];
  TLegend* leg[4];
  const string entry[12] = {
    "x_{A}-x_{B}",
    "x_{B}-x_{C}",
    "x_{C}-x_{A}",
    "y_{A}-y_{B}",
    "y_{B}-y_{C}",
    "y_{C}-y_{A}",
    "x_{A}-x_{B}",
    "x_{B}-x_{C}",
    "x_{C}-x_{A}",
    "y_{A}-y_{B}",
    "y_{B}-y_{C}",
    "y_{C}-y_{A}"
  };

  for (int i=0;i<4;++i) {
    canv[i] = new TCanvas(Form("align_%d",i),"",500,500);
    leg[i] = new TLegend(0.7,0.7,0.9,0.9,"","brNDC");
    leg[i]->SetTextSize(0.05);
    leg[i]->AddEntry(hist[i*3+0],entry[i*3+0].c_str(),"l");
    leg[i]->AddEntry(hist[i*3+1],entry[i*3+1].c_str(),"l");
    leg[i]->AddEntry(hist[i*3+2],entry[i*3+2].c_str(),"l");

    hist[i*3+1]->Draw();
    hist[i*3+2]->Draw("same");
    hist[i*3+0]->Draw("same");
    leg[i]->Draw("same");
    hist[i*3+1]->Draw("axis,same");
    canv[i]->Print(Form("align_%d.png",i));

  }


}

void doMaps(const char*
            dir="/data/users/abelloni/CERN_TB_Aug15_slim_ntuples",
            bool debug=false) {

  // Identify channels we need to use
  struct channel {
    int chan;
    int ieta;
    int idepth;
    string name;
  };

  vector<channel> channels;
  channels.push_back({12,21,4,"SCSN-81"});
  channels.push_back({20,21,6,"EJ-260"});
  channels.push_back({ 9,22,3,"EJ-200_P2"});
  channels.push_back({13,22,4,"EJ-200_2X"});
  channels.push_back({17,22,5,"EJ-200"});

  // Let us define the histogram and book them
  TH2F *hist_eff[5]; // 5 is the number of tiles; channels.size() == 5
  TH2F *hist_den[5];

  TH1F *hist_effX[5]; // 5 is the number of tiles; channels.size() == 5
  TH1F *hist_denX[5];
  TH1F *hist_effY[5]; // 5 is the number of tiles; channels.size() == 5
  TH1F *hist_denY[5];

  if (channels.size()!=5) {
    cout << "Argh, something wrong!" << endl;
    return;
  }
  for (int i=0;i<channels.size();++i) {
    hist_eff[i] = new TH2F(Form("%s_eff",channels[i].name.c_str()),
                           "",350,-75,75,350,-75,75);
    hist_den[i] = new TH2F(Form("%s_den",channels[i].name.c_str()),
                           "",350,-75,75,350,-75,75);
    hist_eff[i]->Sumw2();
    hist_den[i]->Sumw2();

    hist_effX[i] = new TH1F(Form("%s_effX",channels[i].name.c_str()),
                            "",350,-100,100);
    hist_denX[i] = new TH1F(Form("%s_denX",channels[i].name.c_str()),
                            "",350,-100,100);
    hist_effX[i]->Sumw2();
    hist_denX[i]->Sumw2();

    hist_effY[i] = new TH1F(Form("%s_effY",channels[i].name.c_str()),
                            "",400,-100,100);
    hist_denY[i] = new TH1F(Form("%s_denY",channels[i].name.c_str()),
                            "",400,-100,100);
    hist_effY[i]->Sumw2();
    hist_denY[i]->Sumw2();
  }

  // Now we get the data
  TChain* chain = new TChain("slim");
  chain->Add(Form("%s/*_slim.root",dir));

  // Get the branches I need:
  vector<double> *xa = 0,*xb = 0,*xc = 0,*ya = 0,*yb = 0,*yc = 0;
  chain->SetBranchAddress("xA",&xa);
  chain->SetBranchAddress("xB",&xb);
  chain->SetBranchAddress("xC",&xc);
  chain->SetBranchAddress("yA",&ya);
  chain->SetBranchAddress("yB",&yb);
  chain->SetBranchAddress("yC",&yc);

  double pulse[29][10],ped[29];
  chain->SetBranchAddress("pulse",&pulse);
  chain->SetBranchAddress("ped",&ped);

  double intercept_X, slope_X;
  double intercept_Y, slope_Y;
  chain->SetBranchAddress("interceptX",&intercept_X);
  chain->SetBranchAddress("slopeX",&slope_X);
  chain->SetBranchAddress("interceptY",&intercept_Y);
  chain->SetBranchAddress("slopeY",&slope_Y);

  for (unsigned int i=0;i<chain->GetEntries();++i) {
    if (debug && i>1000)
      continue;

    // Do the alignment
    chain->GetEntry(i);

    double x_hit = intercept_X + z_ex*slope_X;
    double y_hit = intercept_Y + z_ex*slope_Y;

    if (i%100000==0)
      std::cout << "  Done with " << i << " / "
                << chain->GetEntriesFast() << " events\n";
    //cout << "(x,y) = (" << x_hit << "," << y_hit << ")\n";

    // loop on the various channels
    for (int i=0;i<channels.size();++i) {

      // Use TS = 5, 6, 7, 8; remove 4 times the pedestal
      double energy_ps =
        pulse[channels[i].chan][5]+
        pulse[channels[i].chan][6]+
        pulse[channels[i].chan][7]+
        pulse[channels[i].chan][8]-
        4*ped[channels[i].chan];

      if (energy_ps>25) {
        hist_eff[i]->Fill(x_hit,y_hit);
        hist_effX[i]->Fill(x_hit);
        hist_effY[i]->Fill(y_hit);
      }
      hist_den[i]->Fill(x_hit,y_hit);
      hist_denX[i]->Fill(x_hit);
      hist_denY[i]->Fill(y_hit);

    } // loop on channels

  } // loop on events

  // here I should plot the efficiency maps, adfter some beautification
  TCanvas *canv[5],*canvX[5],*canvY[5];
  const string entry[5] = {
    "SCSN-81",
    "EJ-260",
    "EJ-200 P2",
    "EJ-200 2X",
    "EJ-200"
  };
  TLine* fid_line = new TLine();
  fid_line->SetLineWidth(2);
  fid_line->SetLineStyle(kDashed);
  for (int i=0;i<channels.size();++i) {
    hist_eff[i]->Divide(hist_eff[i],hist_den[i],1,1,"b");
    hist_eff[i]->GetXaxis()->SetTitle("x [mm]");
    hist_eff[i]->GetYaxis()->SetTitle("y [mm]");

    canv[i] = new TCanvas(channels[i].name.c_str(),"",550,500);
    canv[i]->SetRightMargin(canv[i]->GetLeftMargin());
    hist_eff[i]->Draw("colz");
    fid_line->DrawLine(fiducialX[i][0],fiducialY[i][0],
                       fiducialX[i][1],fiducialY[i][1]);
    fid_line->DrawLine(fiducialX[i][0],fiducialY[i][0],
                       fiducialX[i][2],fiducialY[i][2]);
    fid_line->DrawLine(fiducialX[i][3],fiducialY[i][3],
                       fiducialX[i][1],fiducialY[i][1]);
    fid_line->DrawLine(fiducialX[i][3],fiducialY[i][3],
                       fiducialX[i][2],fiducialY[i][2]);

    TLatex label;
    label.SetNDC();
    label.SetTextSize(0.05);
    label.SetTextAlign(30);
    label.DrawLatex(0.92,0.875,entry[i].c_str());

    canv[i]->Print(Form("efficiency_map_%s.png",channels[i].name.c_str()));
    canv[i]->Print(Form("efficiency_map_%s.C",channels[i].name.c_str()));

    hist_effX[i]->Divide(hist_effX[i],hist_denX[i],1,1,"b");
    hist_effX[i]->GetXaxis()->SetTitle("x [mm]");
    hist_effX[i]->GetYaxis()->SetTitle("Efficiency");

    canvX[i] = new TCanvas(Form("effX_%s",channels[i].name.c_str()),
                           "",500,500);
    hist_effX[i]->Draw();

    TLatex labelX;
    labelX.SetNDC();
    labelX.SetTextSize(0.05);
    labelX.SetTextAlign(30);
    //labelX.DrawLatex(0.92,0.875,entry[i].c_str());
    labelX.DrawLatex(0.8,0.875,entry[i].c_str());

    canvX[i]->Print(Form("efficiency_x_%s.png",channels[i].name.c_str()));

    hist_effY[i]->Divide(hist_effY[i],hist_denY[i],1,1,"b");
    hist_effY[i]->GetXaxis()->SetTitle("y [mm]");
    hist_effY[i]->GetYaxis()->SetTitle("Efficiency");

    canvY[i] = new TCanvas(Form("effY_%s",channels[i].name.c_str()),
                           "",500,500);
    hist_effY[i]->Draw();

    TLatex labelY;
    labelY.SetNDC();
    labelY.SetTextSize(0.05);
    labelY.SetTextAlign(30);
    labelY.DrawLatex(0.92,0.875,entry[i].c_str());

    canvY[i]->Print(Form("efficiency_y_%s.png",channels[i].name.c_str()));

  }

}

const double edges[248] = {
  1.58,   4.73,   7.88,   11.0,   14.2,   17.3,   20.5,   23.6,
  26.8,   29.9,   33.1,   36.2,   39.4,   42.5,   45.7,   48.8,
  53.6,   60.1,   66.6,   73.0,   79.5,   86.0,   92.5,   98.9,
  105,    112,    118,    125,    131,    138,    144,    151,
  157,    164,    170,    177,    186,    199,    212,    225,
  238,    251,    264,    277,    289,    302,    315,    328,
  341,    354,    367,    380,    393,    406,    418,    431,
  444,    464,    490,    516,    542,    568,    594,    620,
  645,    670,    695,    720,    745,
  771,    796,    821,    846,    871,    897,    922,    947,
  960,    1010,   1060,   1120,   1170,   1220,   1270,   1320,
  1370,   1430,   1480,   1530,   1580,   1630,   1690,   1740,
  1790,   1840,   1890,   1940,   2020,   2120,   2230,   2330,
  2430,   2540,   2640,   2740,   2850,   2950,   3050,   3150,
  3260,   3360,   3460,   3570,   3670,   3770,   3880,   3980,
  4080,   4240,   4450,   4650,   4860,   5070,   5280,   5490,
  5680,   5880,   6080,   6280,   6480,
  6680,   6890,   7090,   7290,   7490,   7690,   7890,   8090,
  8400,   8810,   9220,   9630,   10000,  10400,  10900,  11300,
  11700,  12100,  12500,  12900,  13300,  13700,  14100,  14500,
  15000,  15400,  15800,  16200,  16800,  17600,  18400,  19300,
  20100,  20900,  21700,  22500,  23400,  24200,  25000,  25800,
  26600,  27500,  28300,  29100,  29900,  30700,  31600,  32400,
  33200,  34400,  36100,  37700,  39400,  41000,  42700,  44300,
  45900,  47600,  49200,  50800,  52500,
  54100,  55700,  57400,  59000,  60600,  62200,  63900,  65500,
  68000,  71300,  74700,  78000,  81400,  84700,  88000,  91400,
  94700,  98100,  101000, 105000, 108000, 111000, 115000, 118000,
  121000, 125000, 128000, 131000, 137000, 145000, 152000, 160000,
  168000, 176000, 183000, 191000, 199000, 206000, 214000, 222000,
  230000, 237000, 245000, 253000, 261000, 268000, 276000, 284000,
  291000, 302000, 316000, 329000, 343000, 356000, 370000, 384000, 398000};


// Energy; time-slice
void doEnergyTS(const char*
                dir="/data/users/abelloni/CERN_TB_Aug15_slim_ntuples",
                bool debug=false) {


  // Identify channels we need to use
  struct channel {
    int chan;
    int ieta;
    int idepth;
    string name;
  };

  vector<channel> channels;
  channels.push_back({12,21,4,"SCSN-81"});
  channels.push_back({20,21,6,"EJ-260"});
  channels.push_back({ 9,22,3,"EJ-200_P2"});
  channels.push_back({13,22,4,"EJ-200_2X"});
  channels.push_back({17,22,5,"EJ-200"});

  // Let us define the histogram and book them
  TH1F *hist_en[5]; // 5 is the number of tiles; channels.size() == 5
  TH1F *hist_ts[5];
  if (channels.size()!=5) {
    cout << "Argh, something wrong!" << endl;
    return;
  }
  for (int i=0;i<channels.size();++i) {
    hist_en[i] = new TH1F(Form("en_%s",channels[i].name.c_str()),"",247,edges);
    hist_ts[i] = new TH1F(Form("ts_%s",channels[i].name.c_str()),"",
                          10,0.5,10.5);
  }

  // Now we get the data
  TChain* chain = new TChain("slim");
  chain->Add(Form("%s/*_slim.root",dir));

  // Get the branches I need:
  vector<double> *xa = 0,*xb = 0,*xc = 0,*ya = 0,*yb = 0,*yc = 0;
  chain->SetBranchAddress("xA",&xa);
  chain->SetBranchAddress("xB",&xb);
  chain->SetBranchAddress("xC",&xc);
  chain->SetBranchAddress("yA",&ya);
  chain->SetBranchAddress("yB",&yb);
  chain->SetBranchAddress("yC",&yc);

  double pulse[29][10],ped[29];
  chain->SetBranchAddress("pulse",&pulse);
  chain->SetBranchAddress("ped",&ped);

  double intercept_X, slope_X;
  double intercept_Y, slope_Y;
  chain->SetBranchAddress("interceptX",&intercept_X);
  chain->SetBranchAddress("slopeX",&slope_X);
  chain->SetBranchAddress("interceptY",&intercept_Y);
  chain->SetBranchAddress("slopeY",&slope_Y);

  for (unsigned int i=0;i<chain->GetEntries();++i) {
    if (debug && i>1000)
      continue;

    // Do the alignment
    chain->GetEntry(i);

    if (i%100000==0)
      std::cout << "  Done with " << i << " / "
                << chain->GetEntriesFast() << " events\n";

    double x_hit = intercept_X + z_ex*slope_X;
    double y_hit = intercept_Y + z_ex*slope_Y;

    // loop on the various channels
    for (int i=0;i<channels.size();++i) {

      if(!isFiducial(i,x_hit,y_hit))
        continue;

      // Use TS = 5, 6, 7, 8; remove 4 times the pedestal
      double energy_ps =
        pulse[channels[i].chan][5]+
        pulse[channels[i].chan][6]+
        pulse[channels[i].chan][7]+
        pulse[channels[i].chan][8]-
        4*ped[channels[i].chan];

      hist_en[i]->Fill(energy_ps);

      for (int t=0;t<10;++t)
        hist_ts[i]->Fill(t+1,
                         pulse[channels[i].chan][t]-
                         ped[channels[i].chan]);

    } // loop on channels

  } // loop on events

  // Used in both types of plots
  const string entry[5] = {
    "SCSN-81",
    "EJ-260",
    "EJ-200 P2",
    "EJ-200 2X",
    "EJ-200"
  };

  // here I make the plot, after some beautification
  TCanvas* canv[5];
  const int color[5] = {
    kBlack,
    kGreen,
    kBlue,
    kRed,
    kViolet
  };
  const int style[5] = {
    kDotted,
    kSolid,
    kSolid,
    kSolid,
    kDashed
  };
  for (int i=0;i<channels.size();++i) {

    hist_en[i]->GetXaxis()->SetTitle("Charge [fC]");
    hist_en[i]->GetYaxis()->SetTitle("Events");
    hist_en[i]->SetLineWidth(2);
    hist_en[i]->SetLineColor(color[i]);

    canv[i] = new TCanvas(channels[i].name.c_str(),"",500,500);
    canv[i]->SetLogx();
    canv[i]->SetLogy();
    hist_en[i]->Draw("colz");

    TLatex label;
    label.SetNDC();
    label.SetTextSize(0.05);
    label.SetTextAlign(30);
    label.DrawLatex(0.92,0.875,entry[i].c_str());

    label.SetTextAlign(11);
    float eff = hist_en[i]->Integral(hist_en[i]->FindBin(25),
                                     hist_en[i]->GetNbinsX())/
      hist_en[i]->GetEntries();
    //float eff_err = TMath::Sqrt(eff*(1-eff)/hist_en[i]->GetEntries());
    eff*=100;
    //eff_err*=100;
    label.DrawLatex(0.20,0.225,Form("#splitline{#epsilon=%4.1f%%}"
                                    "{Mean=%3.1f#pm%3.1ffC}",
                                    eff,
                                    hist_en[i]->GetMean(),
                                    hist_en[i]->GetMeanError()));

    canv[i]->Print(Form("energy_PS_%s.png",channels[i].name.c_str()));

  }

  TCanvas* canv_ts = new TCanvas("ts","",500,500);

  TLegend* leg = new TLegend(0.2,0.7,0.4,0.9,"","brNDC");
  leg->SetTextSize(0.05);

  for (int i=0;i<channels.size();++i) {

    hist_ts[i]->SetLineWidth(2);
    hist_ts[i]->SetLineColor(color[i]);
    hist_ts[i]->SetLineStyle(style[i]);

    hist_ts[i]->GetXaxis()->SetTitle("Time Slice [25ns]");
    hist_ts[i]->GetXaxis()->SetNdivisions(10);
    hist_ts[i]->GetYaxis()->SetTitle("Charge [fC]");
    hist_ts[i]->GetYaxis()->SetTitleOffset(1.6);

    leg->AddEntry(hist_ts[i],entry[i].c_str(),"l");
  }

  hist_ts[4]->Draw("hist");
  hist_ts[1]->Draw("hist,same");
  hist_ts[2]->Draw("hist,same");
  hist_ts[3]->Draw("hist,same");
  hist_ts[0]->Draw("hist,same");
  leg->Draw("same");
  hist_ts[4]->Draw("axis,same");

  canv_ts->Print("ts.png");

}


// delta-t w.r.t. SCSN-81
void doTime(const char*
            dir="/data/users/abelloni/CERN_TB_Aug15_slim_ntuples",
            bool debug=false) {

  // Identify channels we need to use
  struct channel {
    int chan;
    int ieta;
    int idepth;
    string name;
  };

  vector<channel> channels;
  channels.push_back({12,21,4,"SCSN-81"});
  channels.push_back({20,21,6,"EJ-260"});
  channels.push_back({ 9,22,3,"EJ-200_P2"});
  channels.push_back({13,22,4,"EJ-200_2X"});
  channels.push_back({17,22,5,"EJ-200"});

  // Let us define the histogram and book them
  TH1F *hist_dt[5]; // 5 is the number of tiles; channels.size() == 5
  if (channels.size()!=5) {
    cout << "Argh, something wrong!" << endl;
    return;
  }
  for (int i=0;i<channels.size();++i) {
    hist_dt[i] = new TH1F(Form("dt_%s",channels[i].name.c_str()),"",
                          320,-160,160);
  }

  // Now we get the data
  TChain* chain = new TChain("slim");
  chain->Add(Form("%s/*_slim.root",dir));

  // Get the branches I need:
  vector<double> *xa = 0,*xb = 0,*xc = 0,*ya = 0,*yb = 0,*yc = 0;
  chain->SetBranchAddress("xA",&xa);
  chain->SetBranchAddress("xB",&xb);
  chain->SetBranchAddress("xC",&xc);
  chain->SetBranchAddress("yA",&ya);
  chain->SetBranchAddress("yB",&yb);
  chain->SetBranchAddress("yC",&yc);

  double pulse_tdc[29][10];
  chain->SetBranchAddress("pulse_tdc",&pulse_tdc);

  double pulse[29][10],ped[29];
  chain->SetBranchAddress("pulse",&pulse);
  chain->SetBranchAddress("ped",&ped);

  double intercept_X, slope_X;
  double intercept_Y, slope_Y;
  chain->SetBranchAddress("interceptX",&intercept_X);
  chain->SetBranchAddress("slopeX",&slope_X);
  chain->SetBranchAddress("interceptY",&intercept_Y);
  chain->SetBranchAddress("slopeY",&slope_Y);

  for (unsigned int i=0;i<chain->GetEntries();++i) {
    if (debug && i>1000)
      continue;

    // Do the alignment
    chain->GetEntry(i);

    if (i%100000==0)
      std::cout << "  Done with " << i << " / "
                << chain->GetEntriesFast() << " events\n";

    double x_hit = intercept_X + z_ex*slope_X;
    double y_hit = intercept_Y + z_ex*slope_Y;

    // time_ref is the SCSN-81 time, used in all delta(t)
    double time_ref = -99;

    // loop on the various channels
    for (int i=0;i<channels.size();++i) {

      if(!isFiducial(i,x_hit,y_hit))
        continue;

      // Use TS = 5, 6, 7, 8; remove 4 times the pedestal
      double energy_ps =
	pulse[channels[i].chan][5]+
	pulse[channels[i].chan][6]+
	pulse[channels[i].chan][7]+
	pulse[channels[i].chan][8]-
	4*ped[channels[i].chan];
      if (energy_ps<=25) // require >25fC hits in both scintillators
	continue;
      
      double time = -99;
      // Get best time for current channel
      for (int ts=0;ts<10;++ts)
        if (pulse_tdc[channels[i].chan][ts]>-99 &&
            fabs(pulse_tdc[channels[i].chan][ts])>0.01) {
          time = pulse_tdc[channels[i].chan][ts];
          if(i==0)
            time_ref = time;
        }

      // Now I have the good time, can plot difference
      // Plot will naturally be 0 in SCSN case)
      if (time>-99 && time_ref>-99)
        hist_dt[i]->Fill(time-time_ref);

    } // loop on channels

  } // loop on events

  // here I make the plot, after some beautification
  TCanvas* canv[5];
  const int color[5] = {
    kBlack,
    kGreen,
    kBlue,
    kRed,
    kViolet
  };
  const string entry[5] = {
    "SCSN-81",
    "EJ-260",
    "EJ-200 P2",
    "EJ-200 2X",
    "EJ-200"
  };

  for (int i=0;i<channels.size();++i) {

    hist_dt[i]->GetXaxis()->
      SetTitle(Form("t_{%s}-t_{SCSN-81} [ns]",entry[i].c_str()));
    hist_dt[i]->GetYaxis()->SetTitle("Events");
    hist_dt[i]->SetLineWidth(2);
    hist_dt[i]->SetLineColor(color[i]);

    hist_dt[i]->SetMaximum(hist_dt[i]->GetMaximum()*
                           TMath::Power(10,
                                        TMath::Log(hist_dt[i]->GetMaximum())/
                                        TMath::Log(10)/10.*2));

    canv[i] = new TCanvas(channels[i].name.c_str(),"",500,500);
    canv[i]->SetLogy();
    hist_dt[i]->Draw();

    TLatex label;
    label.SetNDC();
    label.SetTextSize(0.05);
    label.DrawLatex(0.20,0.875,Form("<#deltat>=%3.1f#pm%3.1fns",
                                    hist_dt[i]->GetMean(),
                                    hist_dt[i]->GetMeanError()));

    canv[i]->Print(Form("time_%s.png",channels[i].name.c_str()));

  }


}

bool isFiducial(int i, float x_hit, float y_hit) {

  // Find if it is fiducial!
  int polyCorners = 4;
  int k, j=polyCorners-1 ;
  bool oddNodes=false;
  for (k=0; k<polyCorners; k++) {
    if (((fiducialY[i][k]< y_hit && fiducialY[i][j]>=y_hit) ||
         (fiducialY[i][j]< y_hit && fiducialY[i][k]>=y_hit)) &&
        (fiducialX[i][k]<=x_hit || fiducialX[i][j]<=x_hit)) {
      oddNodes^=(fiducialX[i][k]+(y_hit-fiducialY[i][k])/
                 (fiducialY[i][j]-fiducialY[i][k])*
                 (fiducialX[i][j]-fiducialX[i][k])<x_hit);
    }
    j=k;
  }

  return oddNodes;

}
