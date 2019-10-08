#include "sliman2017.h"

////////////////////////////////////////////////////////////////////////////////
// Alignment plots
////////////////////////////////////////////////////////////////////////////////
void doAlignmentPlots(bool debug, const char* dir) {

  // The idea is to use this function to find the mean position
  // of the wire chambers; that mean position shall then be
  // hard-coded in the slimmer macro, before re-running
  // the slimming jobs and obtain the best track slope

  TChain* chain = new TChain("slim");
  chain->Add(Form("%s/*_slim.root",dir));

  cout << "entries: " << chain->GetEntries() << endl;

  // Plots without correction
  const int NPLANES=4; // uhm, technically it is something different...
  TH1F* hist[NPLANES]; // There will be four histograms
  const string var[NPLANES] = {
    "xC-xA",
    "yC-yA",
    Form("xC-xA%s%f",wc_xca>0?"+":"",wc_xca),
    Form("yC-yA%s%f",wc_yca>0?"+":"",wc_yca),
  };
  const int color[NPLANES] = {
    kBlue,
    kBlue,
    kBlue,
    kBlue
  };

  double mean[NPLANES], sigma[NPLANES];
  for (int i = 0; i < NPLANES; ++i) {
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
    if (var[i].find('x')!=std::string::npos)
      hist[i]->GetXaxis()->SetTitle("#deltax [mm]");
    else
      hist[i]->GetXaxis()->SetTitle("#deltay [mm]");
    hist[i]->GetXaxis()->SetTitleOffset(1.2);
    hist[i]->GetYaxis()->SetTitle("Events");
    hist[i]->GetYaxis()->SetTitleOffset(1.6);
  }

  for (int i = 0; i < NPLANES; ++i)
    cout << var[i].c_str() << " mean  -> " << mean[i] << endl;
  for (int i = 0; i < NPLANES; ++i)
    cout << var[i].c_str() << " sigma -> " << sigma[i] << endl;

  TCanvas* canv[NPLANES];
  TLegend* leg[NPLANES];
  const string entry[NPLANES] = {
    "x_{C}-x_{A}",
    "y_{C}-y_{A}",
    "x_{C}-x_{A}",
    "y_{C}-y_{A}"
  };

  for (int i = 0; i < NPLANES; ++i) {
    canv[i] = new TCanvas(Form("align_%d",i),"",500,500); // makes the new TCanvases                                                            
    leg[i] = new TLegend(0.7,0.7,0.9,0.9,"","brNDC");
    leg[i]->SetTextSize(0.05);
    leg[i]->AddEntry(hist[i],entry[i].c_str(),"l");

    hist[i]->Draw();
    leg[i]->Draw("same");
    hist[i]->Draw("axis,same");
    canv[i]->Print(Form("Alignment_Plots/align_%d.png",i));
    canv[i]->Print(Form("Alignment_Plots/align_%d.pdf",i));
  }
}

// debug = false as default
void doMaps(bool debug, const char* dir) {
  // Identify channels we need to use
  struct channel {
    int chan;
    int ieta;
    int idepth;
    string name;
  };

  bool eff, eff_rot, effX, effY, effX_rot, effY_rot, effX_cut, effY_cut, 
    effX_rot_cut, effY_rot_cut, effX_rot_cut_nbins, effY_rot_cut_nbins, cmb;
  bool crudtest, denXY_rot_cut, pedXY_rot_cut, noisetest, overlay, center;
  debug = true;
  effX_cut = effY_cut = !debug;
  effX_rot_cut = effY_rot_cut = !debug;
  effX_rot_cut_nbins = effY_rot_cut_nbins = !debug;
  effX = effY = !debug;
  effX_rot = effY_rot = !debug;
  // Cmb needs rot to work
  eff_rot = cmb = debug;
  eff = !debug;
  denXY_rot_cut = !debug;
  overlay = false; // true;
  crudtest = false;
  noisetest = false;
  center = false;
  pedXY_rot_cut = !debug;

  // Fill the Rotated Arrays
  // This fills both the Theta array, and the fiducial arrays
  fill_Rot_Array();

  // Let us define the histogram and book them
  TH2F *hist_eff[NUMCHAN];
  TH2F *hist_den[NUMCHAN];
  TH2F *hist_eff_cmb[NUMCHAN];
  TH1F *hist_effX[NUMCHAN]; 
  TH1F *hist_denX[NUMCHAN];
  TH1F *hist_effX_cut[NUMCHAN];
  TH1F *hist_denX_cut[NUMCHAN];
  TH1F *hist_effY[NUMCHAN];
  TH1F *hist_denY[NUMCHAN];
  TH1F *hist_effY_cut[NUMCHAN];
  TH1F *hist_denY_cut[NUMCHAN];

  // Rotated Maps
  TH2F *hist_eff_rot[NUMCHAN]; 
  TH2F *hist_den_rot[NUMCHAN];
  TH1F *hist_effX_rot[NUMCHAN];
  TH1F *hist_denX_rot[NUMCHAN];
  TH1F *hist_effX_rot_cut[NUMCHAN];
  TH1F *hist_denX_rot_cut[NUMCHAN];
  TH1F *hist_effY_rot[NUMCHAN];
  TH1F *hist_denY_rot[NUMCHAN];
  TH1F *hist_effY_rot_cut[NUMCHAN];
  TH1F *hist_denY_rot_cut[NUMCHAN];

  // Special binning
  int bins = 100;
  TH1F *hist_effX_rot_cut_nbins[NUMCHAN];
  TH1F *hist_denX_rot_cut_nbins[NUMCHAN];
  TH1F *hist_effY_rot_cut_nbins[NUMCHAN];
  TH1F *hist_denY_rot_cut_nbins[NUMCHAN];
  
  TH1F *hist_pedX_rot_cut[NUMCHAN];
  TH1F *hist_pedY_rot_cut[NUMCHAN];

  // TEST
  TH2F *hist_testX[NUMCHAN];
  TH2F *hist_testY[NUMCHAN];
  TH1F *hist_noiseY[NUMCHAN];

  
  if (channels.size()!= NUMCHAN) {
    cout << "Argh, something wrong!" << endl;
    return;
  }

  for (unsigned int i = 0; i < channels.size(); ++i) {
    // Horrible way to replace - with _ in the channel's name, so that we can print canvas to .C file
    // and obtain a usable macro
    hist_eff[i] = new TH2F(TString(Form("%s_eff",channels[i].name.c_str())).ReplaceAll("-","_").Data(),
                           "",350,-75,75,350,-75,75);
    hist_eff_cmb[i] = new TH2F(TString(Form("%s_eff_cmb",channels[i].name.c_str())).ReplaceAll("-","_").Data(),
                           "",350,-75,75,350,-75,75);
    hist_den[i] = new TH2F(TString(Form("%s_den",channels[i].name.c_str())).ReplaceAll("-","_").Data(),
                           "",350,-75,75,350,-75,75);
    hist_eff[i]->Sumw2();
    hist_den[i]->Sumw2();

    hist_effX[i] = new TH1F(TString(Form("%s_effX",channels[i].name.c_str())).ReplaceAll("-","_").Data(),
                            "",350, -60, 40);
    hist_denX[i] = new TH1F(TString(Form("%s_denX",channels[i].name.c_str())).ReplaceAll("-","_").Data(),
                            "",350, -60, 40);
    hist_effX_cut[i] = new TH1F(TString(Form("%s_effX_cut",channels[i].name.c_str())).ReplaceAll("-","_").Data(),
                            "",350, -60, 40);
    hist_denX_cut[i] = new TH1F(TString(Form("%s_denX_cut",channels[i].name.c_str())).ReplaceAll("-","_").Data(),
                            "",350, -60, 40);
    hist_effX[i]->Sumw2();
    hist_denX[i]->Sumw2();

    hist_effY[i] = new TH1F(TString(Form("%s_effY",channels[i].name.c_str())).ReplaceAll("-","_").Data(),
                            "",400,-100,100);
    hist_denY[i] = new TH1F(TString(Form("%s_denY",channels[i].name.c_str())).ReplaceAll("-","_").Data(),
                            "",400,-100,100);
    hist_effY_cut[i] = new TH1F(TString(Form("%s_effY_cut",channels[i].name.c_str())).ReplaceAll("-","_").Data(),
                            "",400,-100,100);
    hist_denY_cut[i] = new TH1F(TString(Form("%s_denY_cut",channels[i].name.c_str())).ReplaceAll("-","_").Data(),
                            "",400,-100,100);
    hist_effY[i]->Sumw2();
    hist_denY[i]->Sumw2();

    // Rotated Plots
    // 2D 
    hist_eff_rot[i] = new TH2F(TString(Form("%s_eff_rot",channels[i].name.c_str())).ReplaceAll("-","_").Data(),
                           "",350,-75,75,350,-75,75);
    hist_den_rot[i] = new TH2F(TString(Form("%s_den_rot",channels[i].name.c_str())).ReplaceAll("-","_").Data(),
                           "",350,-75,75,350,-75,75);
    hist_eff_rot[i]->Sumw2();
    hist_den_rot[i]->Sumw2();

    // X Efficiency
    hist_effX_rot[i] = new TH1F(TString(Form("%s_effX_rot",channels[i].name.c_str())).ReplaceAll("-","_").Data(),
				"",350, -60, 40);
    hist_denX_rot[i] = new TH1F(TString(Form("%s_denX_rot",channels[i].name.c_str())).ReplaceAll("-","_").Data(),
				"",350, -60, 40);
    hist_effX_rot_cut[i] = new TH1F(TString(Form("%s_effX_rot_cut",channels[i].name.c_str())).ReplaceAll("-","_").Data(),
				    "",350, -60, 40);
    hist_denX_rot_cut[i] = new TH1F(TString(Form("%s_denX_rot_cut",channels[i].name.c_str())).ReplaceAll("-","_").Data(),
				    "",350, -60, 40);
    hist_effX_rot_cut_nbins[i] = new TH1F(TString(Form("%s_effX_rot_cut_nbins",channels[i].name.c_str())).ReplaceAll("-","_").Data(), "", bins, -60, 60);
    hist_denX_rot_cut_nbins[i] = new TH1F(TString(Form("%s_denX_rot_cut_nbins",channels[i].name.c_str())).ReplaceAll("-","_").Data(), "", bins, -60, 60);

    hist_effX_rot[i]->Sumw2();
    hist_denX_rot[i]->Sumw2();
    hist_effX_rot_cut_nbins[i]->Sumw2();
    hist_denX_rot_cut_nbins[i]->Sumw2();

    // Y Efficiency
    hist_effY_rot[i] = new TH1F(TString(Form("%s_effY_rot",channels[i].name.c_str())).ReplaceAll("-","_").Data(),
				"",400,-60,50);
    hist_denY_rot[i] = new TH1F(TString(Form("%s_denY_rot",channels[i].name.c_str())).ReplaceAll("-","_").Data(),
				"",400,-60,50);
    hist_effY_rot_cut[i] = new TH1F(TString(Form("%s_effY_rot_cut",channels[i].name.c_str())).ReplaceAll("-","_").Data(),
				    "",400,-60,50);
    hist_denY_rot_cut[i] = new TH1F(TString(Form("%s_denY_rot_cut",channels[i].name.c_str())).ReplaceAll("-","_").Data(),
				    "",400,-60,50);
    hist_effY_rot_cut_nbins[i] = new TH1F(TString(Form("%s_effY_rot_cut_nbins",channels[i].name.c_str())).ReplaceAll("-","_").Data(),"", bins,-70, 50);
    hist_denY_rot_cut_nbins[i] = new TH1F(TString(Form("%s_denY_rot_cut_nbins",channels[i].name.c_str())).ReplaceAll("-","_").Data(),"", bins,-70, 50);
    
    hist_effY_rot[i]->Sumw2();
    hist_denY_rot[i]->Sumw2();
    hist_effY_rot_cut_nbins[i]->Sumw2();
    hist_denY_rot_cut_nbins[i]->Sumw2();
    
    // Test Plots
    hist_testX[i] = new TH2F(TString(Form("%s_XTEST",channels[i].name.c_str())).ReplaceAll("-","_").Data(),
			     "",350,-75,75,350,-75,75);
    hist_testY[i] = new TH2F(TString(Form("%s_YTEST",channels[i].name.c_str())).ReplaceAll("-","_").Data(),
			     "",350,-75,75,350,-75,75);
    if (i > 2 || i < 7) {
      hist_noiseY[i] = new TH1F(TString(Form("%s_noiseY",channels[i].name.c_str())).ReplaceAll("-","_").Data(),
				"", 350, -60, 40);
    }
    // Pedestal Plots
    hist_pedX_rot_cut[i] = new TH1F(TString(Form("%s_pedX_rot_cut",channels[i].name.c_str())).ReplaceAll("-","_").Data(),
				    "", 350, -60, 40);
    hist_pedY_rot_cut[i] = new TH1F(TString(Form("%s_ped_rot",channels[i].name.c_str())).ReplaceAll("-","_").Data(),
				    "", bins, -60, 50);
  }

  // Now we get the data
  TChain* chain = new TChain("slim");
  chain->Add(Form("%s/*_slim.root", dir));

  // Get the branches I need:
  vector<double> *xa = 0, *xc = 0,*ya = 0, *yc = 0;
  chain->SetBranchAddress("xA", &xa);
  chain->SetBranchAddress("xC", &xc);
  chain->SetBranchAddress("yA", &ya);
  chain->SetBranchAddress("yC", &yc);

  // Pulse and ped tell us how much we got total and how much was pedestal?
  // Unclear why we needed to perform a pedestal analysis if this data was
  // available
  double pulse[NCH][NTS], ped[NCH];
  chain->SetBranchAddress("pulse", &pulse);
  chain->SetBranchAddress("ped", &ped);

  double intercept_X, slope_X;
  double intercept_Y, slope_Y;
  chain->SetBranchAddress("interceptX", &intercept_X);
  chain->SetBranchAddress("slopeX", &slope_X);
  chain->SetBranchAddress("interceptY", &intercept_Y);
  chain->SetBranchAddress("slopeY", &slope_Y);

  int run;
  chain->SetBranchAddress("run", &run);

  // Event Loop
  int counter = 0;
  int missed = 0;
  // There are 1720742 events (muons)
  for (unsigned int j = 0; j < chain->GetEntries(); ++j) {
    if (!debug && j > 1000)
      continue;
    
    // Do the alignment
    // This fills the variables that had addresses set previously
    chain->GetEntry(j);
    
    // Get the locations of the hits from the wire chambers
    double x_hit = intercept_X + z_ex*slope_X;
    double y_hit = intercept_Y + z_ex*slope_Y;
    
    if (j%100000==0)
      std::cout << "  Done with " << j << " / "
                << chain->GetEntriesFast() << " events\n";
    //cout << "(x,y) = (" << x_hit << "," << y_hit << ")\n";
    
    // loop on the various channels
    for (unsigned int i = 0; i < channels.size(); ++i) {
      // For Finger tiles
      if (run < 3410 && i >= 3 && i <= 6) {
	continue;
      }
      
      double x_hit_rot = rotate_Point(x_hit, y_hit, i, 'X');
      double y_hit_rot = rotate_Point(x_hit, y_hit, i, 'Y');
      
      if ( isFiducial( i, x_hit, y_hit) ^ isRotFiducial( i, x_hit_rot, y_hit_rot)) {
	cout << "ISSUE WITH" << "\n" 
	     << "i = " << i << "\n"
	     << "x_hit = " << x_hit << "\n"
	     << "y_hit = " << y_hit << "\n"
	     << endl;
	counter++;
      }
      
      double energy_ps = 0;
      for (auto ts : TIMESLICES)
        energy_ps += pulse[channels[i].chan][ts];
      energy_ps -= TIMESLICES.size()*ped[channels[i].chan];
      
      if (energy_ps>25) {
        hist_eff[i]->Fill(x_hit,y_hit);
        hist_effX[i]->Fill(x_hit);
        hist_effY[i]->Fill(y_hit);
	// Rotated
	hist_eff_rot[i]->Fill(x_hit_rot, y_hit_rot);
	hist_effX_rot[i]->Fill(x_hit_rot);
	hist_effY_rot[i]->Fill(y_hit_rot);
	// Cut the X Crud
	if ( above_or_below( x_hit_rot, y_hit_rot, i)) {
	  hist_effX_cut[i]->Fill(x_hit);
	  hist_effX_rot_cut[i]->Fill(x_hit_rot);
	  if (center) {
	    if (i == 3)
	      hist_effX_rot_cut_nbins[i]->Fill(x_hit_rot + 32);
	    else if (i == 4)
	      hist_effX_rot_cut_nbins[i]->Fill(x_hit_rot + 15);
	    else if (i == 6)
	      hist_effX_rot_cut_nbins[i]->Fill(x_hit_rot - 17);
	    else
	      hist_effX_rot_cut_nbins[i]->Fill(x_hit_rot);
	  }
	  else
	    hist_effX_rot_cut_nbins[i]->Fill(x_hit_rot);
	    
	  if (j%10 == 0 && crudtest)
	    hist_testX[i] -> Fill(x_hit_rot, y_hit_rot);
	}
	// Cut the Y Crud
	if ( left_or_right( x_hit_rot, y_hit_rot, i)) {
	  hist_effY_cut[i]->Fill(y_hit);
	  hist_effY_rot_cut[i]->Fill(y_hit_rot);
	  hist_effY_rot_cut_nbins[i]->Fill(y_hit_rot);
	  if (j%10 == 0 && crudtest)
	    hist_testY[i] -> Fill(x_hit_rot, y_hit_rot);
	} 
      } // if (energy_ps > 25)
      else {
	if ( above_or_below( x_hit_rot, y_hit_rot, i))
	  hist_pedX_rot_cut[i]->Fill(x_hit_rot);
	if ( left_or_right( x_hit_rot, y_hit_rot, i))
	  hist_pedY_rot_cut[i]->Fill(y_hit_rot);
	if (isOtherFiducial( i, x_hit_rot, y_hit_rot, anti_fiducialX, anti_fiducialY) && noisetest)
	  hist_noiseY[i]->Fill(y_hit_rot);
      }

      // Regular
      hist_den[i]->Fill(x_hit,y_hit);
      hist_denX[i]->Fill(x_hit);
      hist_denY[i]->Fill(y_hit);
      
      // Rotated
      hist_den_rot[i]->Fill(x_hit_rot, y_hit_rot);
      hist_denX_rot[i]->Fill(x_hit_rot);
      hist_denY_rot[i]->Fill(y_hit_rot);

      // Cut denominator plots
      if ( above_or_below( x_hit_rot, y_hit_rot, i)) {
	hist_denX_rot_cut[i]->Fill(x_hit_rot);
	if (center) {
	  if (i == 3)
	    hist_denX_rot_cut_nbins[i]->Fill(x_hit_rot + 32);
	  else if (i == 4)
	    hist_denX_rot_cut_nbins[i]->Fill(x_hit_rot + 15);
	  else if (i == 6)
	    hist_denX_rot_cut_nbins[i]->Fill(x_hit_rot - 17);
	  else
	    hist_denX_rot_cut_nbins[i]->Fill(x_hit_rot);
	}
	else
	  hist_denX_rot_cut_nbins[i]->Fill(x_hit_rot);
      }
      
      if ( left_or_right( x_hit_rot, y_hit_rot, i)) {
	hist_denY_rot_cut[i]->Fill(y_hit_rot);
	hist_denY_rot_cut_nbins[i]->Fill(y_hit_rot);
      }
    } // loop on channels
  } // loop on events
  
  cout << "counter = " << counter << endl;
  cout << "missed for fingers = " << missed << endl;
  
  // here I should plot the efficiency maps, after some beautification
  TCanvas *canv[NUMCHAN], *canvX[NUMCHAN], *canvY[NUMCHAN];
  // Alternate color palette for the 2D Eff
  TCanvas *canv_cmb[NUMCHAN];
  // Plot the rotated hists on these canvases
  TCanvas *canv_rot[NUMCHAN], *canvX_rot[NUMCHAN], *canvY_rot[NUMCHAN];
  // Test Canvases for seeing if cutting crud works
  TCanvas *canvT_Xrot[NUMCHAN], *canvT_Yrot[NUMCHAN];
  // Test Canvas for nouse
  TCanvas *canv_noiseY[NUMCHAN];
  // For cutting crud
  TCanvas *canvX_cut[NUMCHAN], *canvY_cut[NUMCHAN], *canvX_rot_cut[NUMCHAN], *canvY_rot_cut[NUMCHAN];
  // Special Binned Canvases
  TCanvas *canvX_rot_cut_nbins[NUMCHAN], *canvY_rot_cut_nbins[NUMCHAN];
  // Denominator Canvases
  TCanvas *canvdenX_rot_cut[NUMCHAN], *canvdenY_rot_cut[NUMCHAN];
  // Pedestal Canvases
  TCanvas *canvpedX_rot_cut[NUMCHAN], *canvpedY_rot_cut[NUMCHAN];
  // To overlay the fingers
  TCanvas *canv_allfingersY, *canv_allfingersX, *canv_allsigsY, *canv_allsigsX;
  
  for (unsigned int i = 0; i < channels.size(); ++i) {
    // ******** 2D EFFICIENCY ********
    // Set up line
    TLine* fid_line = new TLine();
    fid_line->SetLineWidth(2);
    fid_line->SetLineStyle(kDashed);
    
    // Set Palette to kBird for UMD Style
    gStyle->SetPalette(kBird);
    
    // Divide to get efficiency
    hist_eff[i]->Divide(hist_eff[i],hist_den[i],1,1,"b");
    // Set Axis
    hist_eff[i]->GetXaxis()->SetTitle("x [mm]");
    hist_eff[i]->GetYaxis()->SetTitle("y [mm]");

    if ( eff) {    
      canv[i] = new TCanvas(TString(channels[i].name.c_str()).ReplaceAll("-","_").Data(), "", 550, 500);
      canv[i]->SetRightMargin(canv[i]->GetLeftMargin());

      hist_eff[i]->Draw("colz");
      // Draw Lines
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
      label.DrawLatex(0.8,0.875, entry[i].c_str());
      gPad->SetGrid();
      
      canv[i]->Print(Form("Original_Images/Efficiency_Maps_2D/efficiency_map_%s.png",channels[i].name.c_str()));
      canv[i]->Print(Form("Original_Images/Efficiency_Maps_2D/efficiency_map_%s.pdf",channels[i].name.c_str()));
      canv[i]->Print(Form("Original_Images/Efficiency_Maps_2D/efficiency_map_%s.C",channels[i].name.c_str()));
    }

    // ******* 2D ROTATED EFFICIENCY *******
    TLine* fid_line_rot = new TLine();
    fid_line_rot -> SetLineWidth(2);
    fid_line_rot -> SetLineStyle(kDashed);
    
    hist_eff_rot[i] -> Divide( hist_eff_rot[i], hist_den_rot[i],1,1,"b");
    hist_eff_rot[i] -> GetXaxis() -> SetTitle("x [mm]");
    hist_eff_rot[i] -> GetYaxis() -> SetTitle("y [mm]");

    if ( eff_rot) {
      canv_rot[i] = new TCanvas( TString((channels[i].name + "_rot").c_str()).ReplaceAll("-","_").Data(), "", 550, 500);
      canv_rot[i] -> SetRightMargin(canv_rot[i] -> GetLeftMargin());

      hist_eff_rot[i] -> Draw("colz");
      fid_line_rot -> DrawLine( rot_fiducialX[i][0], rot_fiducialY[i][0],
				rot_fiducialX[i][1], rot_fiducialY[i][1]);
      fid_line_rot -> DrawLine( rot_fiducialX[i][0], rot_fiducialY[i][0],
				rot_fiducialX[i][2], rot_fiducialY[i][2]);
      fid_line_rot -> DrawLine( rot_fiducialX[i][3], rot_fiducialY[i][3],
				rot_fiducialX[i][1], rot_fiducialY[i][1]);
      fid_line_rot -> DrawLine( rot_fiducialX[i][3], rot_fiducialY[i][3],
				rot_fiducialX[i][2], rot_fiducialY[i][2]);

      /*
      fid_line_rot->DrawLine( anti_fiducialX[i][0], anti_fiducialY[i][0],
			      anti_fiducialX[i][1], anti_fiducialY[i][1]);
      fid_line_rot->DrawLine( anti_fiducialX[i][0], anti_fiducialY[i][0],
			      anti_fiducialX[i][2], anti_fiducialY[i][2]);
      fid_line_rot->DrawLine( anti_fiducialX[i][3], anti_fiducialY[i][3],
			      anti_fiducialX[i][1], anti_fiducialY[i][1]);
      fid_line_rot->DrawLine( anti_fiducialX[i][3], anti_fiducialY[i][3],
			      anti_fiducialX[i][2], anti_fiducialY[i][2]);
      */

      TLatex label_rot;
      label_rot.SetNDC();
      label_rot.SetTextSize(0.05);
      label_rot.SetTextAlign(30);
      label_rot.DrawLatex(0.8,0.875, entry[i].c_str());

      canv_rot[i] -> Print(Form("Rotated_Images/Efficiency_Maps_2D/efficiency_map_rot%s.png", channels[i].name.c_str()));
      canv_rot[i] -> Print(Form("Rotated_Images/Efficiency_Maps_2D/efficiency_map_rot%s.pdf", channels[i].name.c_str()));
      canv_rot[i] -> Print(Form("Rotated_Images/Efficiency_Maps_2D/efficiency_map_rot%s.C", channels[i].name.c_str()));
    }

    // ******** CMB-LIKE PLOT OF 2D ROTATED EFFICIENCY ********
    if (cmb) {
      // Define Color Palette
      const Int_t num = 9;
      // These are parameters for "kBird", the palette that we use originaly
      Double_t red[num]   = { 0.2082, 0.0592, 0.0780, 0.0232, 0.1802, 0.5301, 0.8186, 0.9956, 0.9764};
      Double_t green[num] = { 0.1664, 0.3599, 0.5041, 0.6419, 0.7178, 0.7492, 0.7328, 0.7862, 0.9832};
      Double_t blue[num]  = { 0.5293, 0.8684, 0.8385, 0.7914, 0.6425, 0.4662, 0.3499, 0.1968, 0.0539};
      Double_t length[num] = { 0.00, 0.86, 0.88, 0.90, 0.92, 0.94, 0.96, 0.98, 1.0};
      Int_t nb = 255;
      TColor::CreateGradientColorTable( num, length, red, green, blue, nb);
      
      canv_cmb[i] = (TCanvas*)canv_rot[i]->DrawClone();
      canv_cmb[i]->Print(Form("Rotated_Images/Efficiency_Maps_2D/CMB_Plots/efficiency_map_rotcmb%s.png",channels[i].name.c_str()));
      canv_cmb[i]->Print(Form("Rotated_Images/Efficiency_Maps_2D/CMB_Plots/efficiency_map_rotcmb%s.pdf",channels[i].name.c_str()));
      canv_cmb[i]->Print(Form("Rotated_Images/Efficiency_Maps_2D/CMB_Plots/efficiency_map_rotcmb%s.C",channels[i].name.c_str()));
    }

    // ******* X EFFICIENCY *******
    hist_effX[i]->Divide(hist_effX[i],hist_denX[i],1,1,"b");
    hist_effX[i]->GetXaxis()->SetTitle("x [mm]");
    hist_effX[i]->GetYaxis()->SetTitle("Efficiency");
    hist_effX[i]->SetMarkerStyle(7);
    
    if ( effX) {
      canvX[i] = new TCanvas(TString(Form("effX_%s", channels[i].name.c_str())).ReplaceAll("-","_").Data(),
			     "", 500, 500);
      hist_effX[i]->Draw();
      
      TLatex labelX;
      labelX.SetNDC();
      labelX.SetTextSize(0.05);
      labelX.SetTextAlign(30);
      labelX.DrawLatex(0.92,0.875, entry[i].c_str());
      
      canvX[i]->Print(Form("Original_Images/Efficiency_Maps_X/efficiency_X_%s.png", channels[i].name.c_str()));
      canvX[i]->Print(Form("Original_Images/Efficiency_Maps_X/efficiency_X_%s.pdf", channels[i].name.c_str()));
      canvX[i]->Print(Form("Original_Images/Efficiency_Maps_X/efficiency_X_%s.C", channels[i].name.c_str()));
    }

    // ******** X EFFICIENCY CUT ********
    hist_effX_cut[i]->Divide(hist_effX_cut[i],hist_denX_cut[i],1,1,"b");
    hist_effX_cut[i]->GetXaxis()->SetTitle("x [mm]");
    hist_effX_cut[i]->GetYaxis()->SetTitle("Efficiency");
    hist_effX_cut[i]->SetMarkerStyle(7);

    if ( effX_cut) {
      canvX_cut[i] = new TCanvas(TString(Form("effX_cut_%s", channels[i].name.c_str())).ReplaceAll("-","_").Data(),
				 "", 500, 500);
      hist_effX_cut[i]->Draw();
      
      TLatex labelX;
      labelX.SetNDC();
      labelX.SetTextSize(0.05);
      labelX.SetTextAlign(30);
      labelX.DrawLatex(0.92,0.875,(entry[i] + " Cut").c_str() );
      
      canvX_cut[i]->Print(Form("Original_Images/Efficiency_Maps_X/No_Crud/efficiency_X_cut_%s.png", channels[i].name.c_str()));
      canvX_cut[i]->Print(Form("Original_Images/Efficiency_Maps_X/No_Crud/efficiency_X_cut_%s.pdf", channels[i].name.c_str()));
      canvX_cut[i]->Print(Form("Original_Images/Efficiency_Maps_X/No_Crud/efficiency_X_cut_%s.C", channels[i].name.c_str()));
    }
    
    // ******* X EFFICIENCY ROTATED ********
    hist_effX_rot[i]->Divide(hist_effX_rot[i], hist_denX_rot[i], 1, 1, "b");
    hist_effX_rot[i]->GetXaxis()->SetTitle("x [mm]");
    hist_effX_rot[i]->GetYaxis()->SetTitle("Efficiency");
    hist_effX_rot[i]->SetMarkerStyle(7);
    
    if ( effX_rot) {
      canvX_rot[i] = new TCanvas(TString(Form("effX_rot_%s", (channels[i].name + "_rot").c_str())).ReplaceAll("-","_").Data(), "", 500, 500);
      hist_effX_rot[i]->Draw();
      
      TLatex labelX_rot;
      labelX_rot.SetNDC();
      labelX_rot.SetTextSize(0.05);
      labelX_rot.SetTextAlign(30);
      labelX_rot.DrawLatex(0.92, 0.875, (entry[i] + " Rot").c_str());
      
      canvX_rot[i]->Print(Form("Rotated_Images/Efficiency_Maps_X/efficiency_X_rot%s.png", channels[i].name.c_str()));
      canvX_rot[i]->Print(Form("Rotated_Images/Efficiency_Maps_X/efficiency_X_rot%s.pdf", channels[i].name.c_str()));
      canvX_rot[i]->Print(Form("Rotated_Images/Efficiency_Maps_X/efficiency_X_rot%s.C", channels[i].name.c_str()));      
    }

    // ******* X EFFICIENCY CUT ROTATED ********
    hist_effX_rot_cut[i]->Divide(hist_effX_rot_cut[i], hist_denX_rot_cut[i], 1, 1, "b");
    hist_effX_rot_cut[i]->GetXaxis()->SetTitle("x [mm]");
    hist_effX_rot_cut[i]->GetYaxis()->SetTitle("Efficiency");
    hist_effX_rot_cut[i]->SetMarkerStyle(7);

    if ( effX_rot_cut) {
      canvX_rot_cut[i] = new TCanvas(TString(Form("effX_rot_cut_%s", (channels[i].name + "_rot").c_str())).ReplaceAll("-","_").Data(), "", 500, 500);
      
      hist_effX_rot_cut[i]->Draw();
      
      TLatex labelX_rot_cut;
      labelX_rot_cut.SetNDC();
      labelX_rot_cut.SetTextSize(0.04);
      labelX_rot_cut.SetTextAlign(21);
      labelX_rot_cut.DrawLatex(0.92, 0.875, (entry[i] + " Cut & Rot").c_str());
      
      canvX_rot_cut[i]->Print(Form("Rotated_Images/Efficiency_Maps_X/No_Crud/efficiency_X_rot_cut_%s.png", channels[i].name.c_str()));
      canvX_rot_cut[i]->Print(Form("Rotated_Images/Efficiency_Maps_X/No_Crud/efficiency_X_rot_cut_%s.pdf", channels[i].name.c_str()));
      canvX_rot_cut[i]->Print(Form("Rotated_Images/Efficiency_Maps_X/No_Crud/efficiency_X_rot_cut_x%s.C", channels[i].name.c_str()));
    }
    
    // ******* X EFFICIENCY CUT ROTATED SPECIAL BINS  ********
    hist_effX_rot_cut_nbins[i]->Divide(hist_effX_rot_cut_nbins[i], hist_denX_rot_cut_nbins[i], 1, 1, "b");
    hist_effX_rot_cut_nbins[i]->GetXaxis()->SetTitle("x [mm]");
    hist_effX_rot_cut_nbins[i]->GetYaxis()->SetTitle("Efficiency");
    hist_effX_rot_cut_nbins[i]->SetMarkerStyle(7);
    hist_effX_rot_cut_nbins[i]->GetYaxis()->SetRangeUser( 0.5, 1.0);
    hist_effX_rot_cut_nbins[i]->GetYaxis()->SetTitleOffset(1.6);
    
    if ( effX_rot_cut_nbins) {
      // Set Errors
      double bin = 0, denom_bin = 0;
      double error = 0;
      for (int j = 0; j < hist_effX_rot_cut_nbins[i]->GetNbinsX(); j++) {
	bin = hist_effX_rot_cut_nbins[i]->GetBinContent( j);
	denom_bin = hist_denX_rot_cut_nbins[i]->GetBinContent( j);
	error = TMath::Sqrt( bin * ( ( 1 - bin) / denom_bin));
	// std::cout << "Error = " << error << std::endl;
	hist_effX_rot_cut_nbins[i]->SetBinError( j, error);
      }

      std::string my_name = Form("effX_rot_cut_nbins_%s", channels[i].name.c_str());
      canvX_rot_cut_nbins[i] = new TCanvas( my_name.c_str(), "", 500, 500);
      // hist_effX_rot_cut_nbins[i]->Draw( "colz");
      hist_effX_rot_cut_nbins[i]->Draw( );
      
      TLatex labelX_rot_cut_nbins;
      labelX_rot_cut_nbins.SetNDC();
      labelX_rot_cut_nbins.SetTextSize(0.04);
      labelX_rot_cut_nbins.SetTextAlign(22);
      std::string label = "#splitline{" + entry[i] + Form(" Cut & Rot}{w/ %d bins", bins) + "}";
      labelX_rot_cut_nbins.DrawLatex(0.5, 0.3, label.c_str());
      
      canvX_rot_cut_nbins[i]->Print(Form("Rotated_Images/Efficiency_Maps_X/No_Crud/Special_Bins/%s.png", my_name.c_str()));
      canvX_rot_cut_nbins[i]->Print(Form("Rotated_Images/Efficiency_Maps_X/No_Crud/Special_Bins/%s.pdf", my_name.c_str()));
      canvX_rot_cut_nbins[i]->Print(Form("Rotated_Images/Efficiency_Maps_X/No_Crud/Special_Bins/%s.C", my_name.c_str()));
    }
    
    // ******* Y EFFICIENCY *******
    hist_effY[i]->Divide(hist_effY[i],hist_denY[i],1,1,"b");
    hist_effY[i]->GetXaxis()->SetTitle("y [mm]");
    hist_effY[i]->GetYaxis()->SetTitle("Efficiency");
    hist_effY[i]->SetMarkerStyle(7);
    
    if ( effY) {
      canvY[i] = new TCanvas(TString(Form("effY_%s",channels[i].name.c_str())).ReplaceAll("-","_").Data(),
			     "",500,500);
      hist_effY[i]->Draw();
      
      TLatex labelY;
      labelY.SetNDC();
      labelY.SetTextSize(0.05);
      labelY.SetTextAlign(30);
      labelY.DrawLatex(0.92,0.875,entry[i].c_str());
      
      canvY[i]->Print(Form("Original_Images/Efficiency_Maps_Y/efficiency_Y_%s.png", channels[i].name.c_str()));
      canvY[i]->Print(Form("Original_Images/Efficiency_Maps_Y/efficiency_Y_%s.pdf", channels[i].name.c_str()));
      canvY[i]->Print(Form("Original_Images/Efficiency_Maps_Y/efficiency_Y_%s.C", channels[i].name.c_str()));
    }

    // ******** Y EFFICIENCY CUT ********
    hist_effY_cut[i]->Divide(hist_effY_cut[i], hist_denY_cut[i], 1, 1, "b");
    hist_effY_cut[i]->GetXaxis()->SetTitle("y [mm]");
    hist_effY_cut[i]->GetYaxis()->SetTitle("Efficiency");
    hist_effY_cut[i]->SetMarkerStyle(7);

    if ( effY_cut) {
      canvY_cut[i] = new TCanvas(TString(Form("effY_cut_%s", channels[i].name.c_str())).ReplaceAll("-","_").Data(),
				 "", 500, 500);
      hist_effY_cut[i]->Draw();
      
      TLatex labelY;
      labelY.SetNDC();
      labelY.SetTextSize(0.05);
      labelY.SetTextAlign(30);
      labelY.DrawLatex(0.92,0.875, (entry[i] + " Cut").c_str());
      
      canvY_cut[i]->Print(Form("Original_Images/Efficiency_Maps_Y/No_Crud/efficiency_Y_cut_%s.png", channels[i].name.c_str()));
      canvY_cut[i]->Print(Form("Original_Images/Efficiency_Maps_Y/No_Crud/efficiency_Y_cut_%s.pdf", channels[i].name.c_str()));
      canvY_cut[i]->Print(Form("Original_Images/Efficiency_Maps_Y/No_Crud/efficiency_Y_cut_%s.C", channels[i].name.c_str()));
    }
    
    // ******* Y EFFICIENCY ROTATED *******
    hist_effY_rot[i]->Divide(hist_effY_rot[i], hist_denY_rot[i], 1, 1, "b");
    hist_effY_rot[i]->GetXaxis()->SetTitle("y [mm]");
    hist_effY_rot[i]->GetYaxis()->SetTitle("Efficiency");
    hist_effY_rot[i]->SetMarkerStyle(7);
    
    if ( effY_rot) {
      canvY_rot[i] = new TCanvas(TString(Form("effY_%s", (channels[i].name + "_rot").c_str())).ReplaceAll("-","_").Data(),"",500,500);
      hist_effY_rot[i]->Draw();
      
      TLatex labelY_rot;
      labelY_rot.SetNDC();
      labelY_rot.SetTextSize(0.05);
      labelY_rot.SetTextAlign(30);
      labelY_rot.DrawLatex(0.92,0.875, (entry[i] + " Rot").c_str());
      
      canvY_rot[i]->Print(Form("Rotated_Images/Efficiency_Maps_Y/efficiency_Y_rot_%s.png", channels[i].name.c_str()));
      canvY_rot[i]->Print(Form("Rotated_Images/Efficiency_Maps_Y/efficiency_Y_rot_%s.pdf", channels[i].name.c_str()));
      canvY_rot[i]->Print(Form("Rotated_Images/Efficiency_Maps_Y/efficiency_Y_rot_%s.C", channels[i].name.c_str()));
    }

    // ******* Y EFFICIENCY CUT ROTATED *******
    hist_effY_rot_cut[i]->Divide(hist_effY_rot_cut[i], hist_denY_rot_cut[i], 1, 1, "b");
    hist_effY_rot_cut[i]->GetXaxis()->SetTitle("y [mm]");
    hist_effY_rot_cut[i]->GetYaxis()->SetTitle("count");
    // if it's not in the fingers
    if (i < 3 || i == 7)
      hist_effY_rot_cut[i]->GetYaxis()->SetRangeUser( 0.5, 1.0);
    hist_effY_rot_cut[i]->SetMarkerStyle(7);
    
    if ( effY_rot_cut) {
      canvY_rot_cut[i] = new TCanvas(TString(Form("effY_%s", (channels[i].name + "_rot").c_str())).ReplaceAll("-","_").Data(),"",500,500);
      hist_effY_rot_cut[i]->Draw();
      
      TLatex labelY_rot_cut;
      labelY_rot_cut.SetNDC();
      labelY_rot_cut.SetTextSize(0.04);
      labelY_rot_cut.SetTextAlign(21);
      labelY_rot_cut.DrawLatex(0.92,0.875, (entry[i] + " Cut & Rot").c_str());
      
      canvY_rot_cut[i]->Print(Form("Rotated_Images/Efficiency_Maps_Y/No_Crud/efficiency_Y_rot_cut_%s.png", channels[i].name.c_str()));
      canvY_rot_cut[i]->Print(Form("Rotated_Images/Efficiency_Maps_Y/No_Crud/efficiency_Y_rot_cut_%s.pdf", channels[i].name.c_str()));
      canvY_rot_cut[i]->Print(Form("Rotated_Images/Efficiency_Maps_Y/No_Crud/efficiency_Y_rot_cut_%s.C", channels[i].name.c_str()));
    }
    // ******** Y EFFICIENCY CUT ROTATED SPECIAL BINS ********
    hist_effY_rot_cut_nbins[i]->Divide(hist_effY_rot_cut_nbins[i], hist_denY_rot_cut_nbins[i], 1, 1, "b");
    hist_effY_rot_cut_nbins[i]->GetXaxis()->SetTitle("y [mm]");
    hist_effY_rot_cut_nbins[i]->GetYaxis()->SetTitle("Efficiency");
    hist_effY_rot_cut_nbins[i]->SetMarkerStyle(7);
    hist_effY_rot_cut_nbins[i]->GetYaxis()->SetRangeUser( 0.5, 1.0);
    hist_effY_rot_cut_nbins[i]->GetYaxis()->SetTitleOffset(1.6);

    if ( effY_rot_cut_nbins) {
      // Set Errors
      double bin = 0, denom_bin = 0;
      double error = 0;
      for (int j = 0; j < hist_effY_rot_cut_nbins[i]->GetNbinsX(); j++) {
	bin = hist_effY_rot_cut_nbins[i]->GetBinContent( j);
	denom_bin = hist_denY_rot_cut_nbins[i]->GetBinContent( j);
	error = TMath::Sqrt( bin * ( ( 1 - bin) / denom_bin));
	// std::cout << "Error = " << error << std::endl;
	hist_effY_rot_cut_nbins[i]->SetBinError( j, error);
      }

      std::string my_name = Form("effY_rot_cut_nbins_%s", channels[i].name.c_str());
      canvY_rot_cut[i] = new TCanvas( my_name.c_str(),"",500,500);
      hist_effY_rot_cut_nbins[i]->Draw( "colz");
      
      TLatex labelY_rot_cut_nbins;
      labelY_rot_cut_nbins.SetNDC();
      labelY_rot_cut_nbins.SetTextSize(0.04);
      labelY_rot_cut_nbins.SetTextAlign(22);
      std::string label = "#splitline{" + entry[i] + Form(" Cut & Rot}{w/ %d bins", bins) + "}";
      labelY_rot_cut_nbins.DrawLatex(0.5, 0.3, label.c_str());

      canvY_rot_cut[i]->Print(Form("Rotated_Images/Efficiency_Maps_Y/No_Crud/Special_Bins/%s.png", my_name.c_str()));
      canvY_rot_cut[i]->Print(Form("Rotated_Images/Efficiency_Maps_Y/No_Crud/Special_Bins/%s.pdf", my_name.c_str()));
      canvY_rot_cut[i]->Print(Form("Rotated_Images/Efficiency_Maps_Y/No_Crud/Special_Bins/%s.C", my_name.c_str()));
    }

    // ******* PRINT DENOMINATORS *******
    if ( denXY_rot_cut) {
      canvdenX_rot_cut[i] = new TCanvas(TString(Form("denX_%s", ((channels[i].name + "_rot").c_str()))).ReplaceAll("-","_").Data(), "", 500, 500);
      
      hist_denX_rot_cut[i]->Draw();
      TLatex label_denomX;
      label_denomX.SetNDC();
      label_denomX.SetTextSize(0.05);
      label_denomX.SetTextAlign(30);
      label_denomX.DrawLatex(0.92,0.875, (entry[i] + Form("X Denominator rot & cut")).c_str());
      gPad->Modified();
      canvdenX_rot_cut[i]->Print(Form("Denominator_Plots/denom_X_rot_cut_%s.png", channels[i].name.c_str()));
      canvdenX_rot_cut[i]->Print(Form("Denominator_Plots/denom_X_rot_cut_%s.pdf", channels[i].name.c_str()));
      
      canvdenY_rot_cut[i] = new TCanvas(TString(Form("denY_%s", ((channels[i].name + "_rot").c_str()))).ReplaceAll("-","_").Data(), "", 500, 500);
      hist_denY_rot_cut[i]->Draw();
      TLatex label_denomY;
      label_denomY.SetNDC();
      label_denomY.SetTextSize(0.05);
      label_denomY.SetTextAlign(30);
      label_denomY.DrawLatex(0.92,0.875, (entry[i] + Form("Y Denominator rot & cut")).c_str());
      gPad->Modified();
      canvdenY_rot_cut[i]->Print(Form("Denominator_Plots/denom_Y_rot_cut_%s.png", channels[i].name.c_str()));
      canvdenY_rot_cut[i]->Print(Form("Denominator_Plots/denom_Y_rot_cut_%s.pdf", channels[i].name.c_str()));
    }

    // ******* PRINT PEDESTAL *******
    hist_pedX_rot_cut[i]->GetXaxis()->SetTitle("x [mm]");
    hist_pedX_rot_cut[i]->GetYaxis()->SetTitle("count");
    hist_pedY_rot_cut[i]->GetXaxis()->SetTitle("y [mm]");
    hist_pedY_rot_cut[i]->GetYaxis()->SetTitle("count"); 
    if ( pedXY_rot_cut) {
      canvpedX_rot_cut[i] = new TCanvas(TString(Form("pedX_%s", ((channels[i].name + "_rot").c_str()))).ReplaceAll("-","_").Data(), "", 500, 500);
      
      hist_pedX_rot_cut[i]->Draw();
      TLatex label_pedX;
      label_pedX.SetNDC();
      label_pedX.SetTextSize(0.04);
      label_pedX.SetTextAlign(22);
      label_pedX.DrawLatex(0.5,0.3, (entry[i] + Form(" X Ped rot & cut")).c_str());
      gPad->Modified();
      canvpedX_rot_cut[i]->Print(Form("Pedestal_Plots/ped_X_rot_cut_%s.png", channels[i].name.c_str()));
      canvpedX_rot_cut[i]->Print(Form("Pedestal_Plots/ped_X_rot_cut_%s.pdf", channels[i].name.c_str()));
      
      canvpedY_rot_cut[i] = new TCanvas(TString(Form("pedY_%s", ((channels[i].name + "_rot").c_str()))).ReplaceAll("-","_").Data(), "", 500, 500);
      
      hist_pedY_rot_cut[i]->Draw();
      TLatex label_pedY;
      label_pedY.SetNDC();
      label_pedY.SetTextSize(0.04);
      label_pedY.SetTextAlign(22);
      label_pedY.DrawLatex(0.5,0.3, (entry[i] + Form(" Y Ped rot & cut")).c_str());
      gPad->Modified();
      canvpedY_rot_cut[i]->Print(Form("Pedestal_Plots/ped_Y_rot_cut_%s.png", channels[i].name.c_str()));
      canvpedY_rot_cut[i]->Print(Form("Pedestal_Plots/ped_Y_rot_cut_%s.pdf", channels[i].name.c_str()));
    }

  
  } // loop on channels (i)

  if (overlay) {
    // ******** OVERLAYED FINGER TILES X ********
    TLegend* leg_fingX = new TLegend(0.2,0.7,0.4,0.9,"","brNDC");
    canv_allfingersX = new TCanvas("Overlayed_FingersX", "", 500, 500);
    
    for (auto k : fingers) {
      hist_effX_rot_cut_nbins[k]->SetLineColor(color[k - 3]);
      hist_effX_rot_cut_nbins[k]->Draw("same, e1");
      hist_effX_rot_cut_nbins[k]->SetMarkerStyle(20);
      hist_effX_rot_cut_nbins[k]->SetMarkerSize(0.5);
      hist_effX_rot_cut_nbins[k]->SetMarkerColor(color[k - 3]);
      hist_effX_rot_cut_nbins[k]->GetYaxis()->SetRangeUser( 0.4, 0.8);
      leg_fingX->AddEntry(hist_effX_rot_cut_nbins[k], entry[k].c_str(), "l");
    }
    hist_effX_rot_cut_nbins[fingers[0]]->Draw("same, axis");

    TLatex label_fingX;
    // NDC means normalized coordinates
    label_fingX.SetNDC();
    label_fingX.SetTextSize(0.05);
    label_fingX.SetTextAlign(30);
    //label_fingX.DrawLatex(0.92,0.875, "Finger Tiles X Eff");

    leg_fingX->SetTextSize(.03);
    leg_fingX->Draw();
  
    // Now update
    gPad->Update();
    if (center) {
      canv_allfingersX->Print("Overlayed_Plots/Overlayed_Finger_XEff_nbins_centered.png");
      canv_allfingersX->Print("Overlayed_Plots/Overlayed_Finger_XEff_nbins_centered.pdf");
      canv_allfingersX->Print("Overlayed_Plots/Overlayed_Finger_XEff_nbins_centered.C");
    }
    else {
      canv_allfingersX->Print("Overlayed_Plots/Overlayed_Finger_XEff_nbins.png");
      canv_allfingersX->Print("Overlayed_Plots/Overlayed_Finger_XEff_nbins.pdf");
      canv_allfingersX->Print("Overlayed_Plots/Overlayed_Finger_XEff_nbins.C");
    }

    // ******** OVERLAYED FINGER TILES Y ********
    TLegend* leg_fingY = new TLegend(0.4,0.2,0.6,0.5,"","brNDC");
    canv_allfingersY = new TCanvas("Overlayed_Fingers", "", 500, 500);
    
    for (auto k : fingers) {
      hist_effY_rot_cut_nbins[k]->SetLineColor(color[k - 3]);
      hist_effY_rot_cut_nbins[k]->Draw("same, e1");
      hist_effY_rot_cut_nbins[k]->SetMarkerStyle(20);
      hist_effY_rot_cut_nbins[k]->SetMarkerSize(0.5);
      hist_effY_rot_cut_nbins[k]->SetMarkerColor(color[k - 3]);
      hist_effY_rot_cut_nbins[k]->GetYaxis()->SetRangeUser( 0.4, 0.75);
      leg_fingY->AddEntry(hist_effY_rot_cut_nbins[k], entry[k].c_str(), "l");
    }
    hist_effY_rot_cut_nbins[fingers[0]]->Draw("same, axis");
    
    TLatex label_fing;
    // NDC means normalized coordinates
    label_fing.SetNDC();
    label_fing.SetTextSize(0.05);
    label_fing.SetTextAlign(30);
    //label_fing.DrawLatex(0.92,0.875, "Finger Tiles Y Eff");
  
    leg_fingY->SetTextSize(.03);
    leg_fingY->Draw();
    // Have to do this after a Draw
  
    // Now update
    gPad->Update();
    canv_allfingersY->Print("Overlayed_Plots/Overlayed_Finger_YEff_nbins.png");
    canv_allfingersY->Print("Overlayed_Plots/Overlayed_Finger_YEff_nbins.pdf");
    canv_allfingersY->Print("Overlayed_Plots/Overlayed_Finger_YEff_nbins.C");

    // ******* OVERLAYED SIGMAS Y BINS ********
    TLegend *leg_sigY = new TLegend(0.4, 0.2, 0.7, 0.4,"","brNDC");
    canv_allsigsY = new TCanvas("Overlayed_SigmasY", "", 500, 500);

    for (auto k : sigmas) {
      hist_effY_rot_cut_nbins[k]->SetLineColor(color[k]);
      hist_effY_rot_cut_nbins[k]->Draw("same, e1");
      hist_effY_rot_cut_nbins[k]->SetMarkerStyle(20);
      hist_effY_rot_cut_nbins[k]->SetMarkerSize(0.5);
      hist_effY_rot_cut_nbins[k]->SetMarkerColor(color[k]);
      hist_effY_rot_cut_nbins[k]->GetYaxis()->SetRangeUser( 0.4, 0.75);
      hist_effY_rot_cut_nbins[k]->GetYaxis()->SetRangeUser( 0.6, 0.9);
      leg_sigY->AddEntry(hist_effY_rot_cut_nbins[k], entry[k].c_str(), "l");
    }
    hist_effY_rot_cut_nbins[sigmas[0]]->Draw("same, axis");

    leg_sigY->Draw();

    TLatex label_sigY;
    // NDC means normalized coordinates
    label_sigY.SetNDC();
    label_sigY.SetTextSize(0.05);
    label_sigY.SetTextAlign(31);
    //label_sigY.DrawLatex(0.92,0.875, "Sigma Tiles Y Eff");
  
    // Now update
    gPad->Update();
    canv_allsigsY->Print("Overlayed_Plots/Overlayed_Sigmas_YEff_nbins.png");
    canv_allsigsY->Print("Overlayed_Plots/Overlayed_Sigmas_YEff_nbins.pdf");
    canv_allsigsY->Print("Overlayed_Plots/Overlayed_Sigmas_YEff_nbins.C");
    
    // ******* OVERLAYED SIGMAS X BINS *******
    TLegend *leg_sigX = new TLegend(0.4, 0.2, 0.7, 0.5,"","brNDC");
    canv_allsigsX = new TCanvas("Overlayed_SigmasX", "", 500, 500);
    
    for (auto k : sigmas) {
      hist_effX_rot_cut_nbins[k]->SetLineColor(color[k]);
      hist_effX_rot_cut_nbins[k]->Draw("same, e1");
      hist_effX_rot_cut_nbins[k]->SetMarkerStyle(20);
      hist_effX_rot_cut_nbins[k]->SetMarkerSize(0.5);
      hist_effX_rot_cut_nbins[k]->SetMarkerColor(color[k]);
      hist_effX_rot_cut_nbins[k]->GetYaxis()->SetRangeUser( 0.4, 0.75);
      hist_effX_rot_cut_nbins[k]->GetYaxis()->SetRangeUser( 0.6, 0.9);
      leg_sigX->AddEntry(hist_effX_rot_cut_nbins[k], entry[k].c_str(), "l");
    }
    hist_effX_rot_cut_nbins[sigmas[0]]->Draw("same, axis");

    leg_sigX->Draw();
    TLatex label_sigX;
    label_sigX.SetNDC();
    label_sigX.SetTextSize(0.05);
    label_sigX.SetTextAlign(31);
    //label_sigX.DrawLatex(0.92,0.875,"Sigma Tiles X Eff");
  
    // Now update
    gPad->Update();
    canv_allsigsX->Print("Overlayed_Plots/Overlayed_Sigmas_XEff_nbins.png");
    canv_allsigsX->Print("Overlayed_Plots/Overlayed_Sigmas_XEff_nbins.pdf");
    canv_allsigsX->Print("Overlayed_Plots/Overlayed_Sigmas_XEff_nbins.C");
  }
  
  // ******** NOISE TEST ********
  if (noisetest) {
    for (unsigned int i = 0; i < NUMCHAN; i++) {
      hist_noiseY[i]->GetXaxis()->SetTitle("y [mm]");
      hist_noiseY[i]->GetYaxis()->SetTitle("count");
      canv_noiseY[i] = new TCanvas(TString(Form("noiseY_%s", channels[i].name.c_str())).ReplaceAll("-","_").Data(),
				   "",500,500);
      hist_noiseY[i]->Draw();

      TLatex label_noise;
      label_noise.SetNDC();
      label_noise.SetTextSize(0.05);
      label_noise.SetTextAlign(31);
      label_noise.DrawLatex(0.92,0.875, (entry[i] + " Noise").c_str());

      // Now update
      gPad->Update();
      canv_noiseY[i]->Print(Form( "Noise_Plots/noiseY_%s.png", channels[i].name.c_str()));
      canv_noiseY[i]->Print(Form( "Noise_Plots/noiseY_%s.pdf", channels[i].name.c_str()));
      delete canv_noiseY[i];
    }
  }
  // ******** CRUD TEST ********
  if (crudtest) {
    for (unsigned int i = 0; i < NUMCHAN; i++) {
      canvT_Xrot[i] = new TCanvas(TString(Form("XTEST_%s", channels[i].name.c_str())).ReplaceAll("-","_").Data(),
				  "",500,500);
      hist_testX[i]->SetMarkerSize(.25);
      hist_testX[i]->GetXaxis()->SetTitle("x [mm]");
      hist_testX[i]->GetYaxis()->SetTitle("y [mm]");
      hist_testX[i]->Draw();

      TLatex label_testX;
      label_testX.SetNDC();
      label_testX.SetTextSize(0.05);
      label_testX.SetTextAlign(31);
      label_testX.DrawLatex( 0.92, 0.875, entry[i].c_str());

      canvT_Xrot[i] -> Print( Form( "Crud_test/crud_Xtest%s.png", channels[i].name.c_str()));
      canvT_Xrot[i] -> Print( Form( "Crud_test/crud_Xtest%s.pdf", channels[i].name.c_str()));
      
      canvT_Yrot[i] = new TCanvas(TString(Form("YTEST_%s", channels[i].name.c_str())).ReplaceAll("-","_").Data(),
				  "",500,500);
      hist_testY[i]->SetMarkerSize(.25);
      hist_testY[i]->GetXaxis()->SetTitle("x [mm]");
      hist_testY[i]->GetYaxis()->SetTitle("y [mm]");
      hist_testY[i]->Draw();
      
      TLatex label_testY;
      label_testY.SetNDC();
      label_testY.SetTextSize(0.05);
      label_testY.SetTextAlign(31);
      label_testY.DrawLatex( 0.92, 0.875, entry[i].c_str());
      
      canvT_Yrot[i] -> Print( Form( "Crud_test/crud_Ytest%s.png", channels[i].name.c_str()));
      canvT_Yrot[i] -> Print( Form( "Crud_test/crud_Ytest%s.pdf", channels[i].name.c_str()));
    }
  }
  
  // If you want to change between regular and other fiducial areas, change the variable in "slimanalysis.h"
#if REG
  cout << "using REGULAR fiducial area" << endl;
#else
  cout << "using OTHER fiducial area" << endl;
#endif
  
} // void doMaps(...


////////////////////////////////////////////////////////////////////////////////
// Energy
////////////////////////////////////////////////////////////////////////////////
void doEnergy(int flag, bool debug, const char* dir) {

  // Make a TFile: will save histograms here, for possible fits later
  TFile *energy_hists = new TFile("energy_hists.root", "RECREATE");
  if (!energy_hists->IsOpen()) {
    cout << "Something went wrong..." << endl;
    return;
  }

  // Default: do all!
  bool reg = 1 & flag;
  bool bins = 2 & flag;
  bool pedestal = 4 & flag;
  bool overlay = 8 & flag;
  
  // Let us define the histogram and book them
  TH1F *hist_en[NUMCHAN]; // 8 is the number of tiles; channels.size() == 8
  TH1F *hist_en_bins[NUMCHAN]; // This will be used for pulse shape fits
  TH1F *hist_en_ped[NUMCHAN]; // This uses a region NOT where the tile is

  if (channels.size()!= NUMCHAN) {
    cout << "Argh, something wrong!" << endl;
    return;
  }

  const unsigned int enbins = 300;
  for (unsigned int i = 0; i < channels.size(); ++i) {
    hist_en[i] = new TH1F( Form("en_%s", channels[i].name.c_str()), "", 247, edges);
    hist_en_bins[i] = new TH1F( Form("en_bins_%s", channels[i].name.c_str()), "", enbins, 0.0, 600.0);
    hist_en_ped[i] = new TH1F( Form("en_ped_%s", channels[i].name.c_str()), "", 37, -50,50); // binning here very fine tuned...
  }

  // Now we get the data
  TChain* chain = new TChain("slim");
  chain->Add(Form("%s/*_slim.root", dir));

  // Get the branches I need:
  vector<double> *xa = 0, *xc = 0, *ya = 0, *yc = 0;
  chain->SetBranchAddress("xA",&xa);
  chain->SetBranchAddress("xC",&xc);
  chain->SetBranchAddress("yA",&ya);
  chain->SetBranchAddress("yC",&yc);

  // Alberto: this was double in James's macro
  float pulse[NCH][NTS], ped[NCH];
  chain->SetBranchAddress("pulse",&pulse);
  chain->SetBranchAddress("ped",&ped);

  double intercept_X, slope_X;
  double intercept_Y, slope_Y;
  chain->SetBranchAddress("interceptX",&intercept_X);
  chain->SetBranchAddress("slopeX",&slope_X);
  chain->SetBranchAddress("interceptY",&intercept_Y);
  chain->SetBranchAddress("slopeY",&slope_Y);

  int run;
  chain->SetBranchAddress("run", &run);

  // Here we fill all histograms
  for (unsigned int i=0;i<chain->GetEntries();++i) {

    // If debugging, let us run only on 50k events
    if (debug && i>50000)
      continue;

    // Do the alignment
    chain->GetEntry(i);

    if (i%100000==0)
      std::cout << "  Done with " << i << " / "
                << chain->GetEntriesFast() << " events\n";

    double x_hit = intercept_X + z_ex*slope_X;
    double y_hit = intercept_Y + z_ex*slope_Y;

    // loop on the various channels
    for (unsigned int i = 0; i < channels.size(); ++i) {
      // This ensures that the finger tiles are using the correct runs
      if (run < 3410 && i >= 3 && i <= 6)
        continue;

      // Use the time slices numbers contained in TIMESLICES
      // and remove pedestal times # of time slices used 
      double energy_ps = 0;
      for (auto ts : TIMESLICES)
        energy_ps+=pulse[channels[i].chan][ts];

      // I think this is the pedestal cut?
      energy_ps-=TIMESLICES.size()*ped[channels[i].chan];

      // Look for anti-fiduciality to fill pedestal plot
      // not enough to say !isFiducial, want to leave some space between fiducial region
      // and anti-fiducial region
      if ( isOtherFiducial( i, x_hit, y_hit, anti_fiducialX, anti_fiducialY))
	hist_en_ped[i]->Fill(energy_ps);
      
      // THIS IS THE MOST IMPORTANT STEP 
      if (!isFiducial( i, x_hit, y_hit))
        continue;

      hist_en[i]->Fill(energy_ps);
      hist_en_bins[i]->Fill(energy_ps);

    } // loop on channels
  } // loop on events

  // All histograms are now filled, we need to make them look good

  // Make energy plots
  TCanvas* canv[NUMCHAN];
  TCanvas* canv_bins[NUMCHAN];
  TCanvas* canv_ped[NUMCHAN];

  for (unsigned int i = 0; i < channels.size(); ++i) {

    hist_en[i]->GetXaxis()->SetTitle("Charge [fC]");
    hist_en[i]->GetYaxis()->SetTitle("Events");
    hist_en[i]->SetLineWidth(2);
    hist_en[i]->SetLineColor(color[i]);
    cout << "After fiducial cut, we use " << (int)hist_en[i]->GetEntries() 
	 << " entries for " << entry[i].c_str() << endl;

    hist_en_ped[i]->GetXaxis()->SetTitle("Charge [fC]");
    hist_en_ped[i]->GetYaxis()->SetTitle("Events");
    hist_en_ped[i]->SetLineWidth(2);
    hist_en_ped[i]->SetLineColor(color[i]);
    cout << "After (anti-)fiducial cut, we use " << (int)hist_en_ped[i]->GetEntries() 
	 << " entries for " << entry[i].c_str() << " pedestal" << endl;

    hist_en_bins[i]->GetXaxis()->SetTitle("Charge [fC]");
    hist_en_bins[i]->GetYaxis()->SetTitle("Events");
    hist_en_bins[i]->SetLineWidth(2);
    hist_en_bins[i]->SetLineColor(color[i]);
    if (hist_en[i]->GetEntries()!=hist_en_bins[i]->GetEntries()) {
      cout << "Argh, the two energy plots should have the same number of entries!" << endl
	   << "---------> Problem in " << entry[i].c_str() << endl;
      return;
    }

    // ********************************************************************************
    // ******** OPTIMIZED BINNING ENERGY HISTOGRAM *******
    // ********************************************************************************
    if (reg) {
      canv[i] = new TCanvas(Form("en_canv_%s", channels[i].name.c_str()), "", CANVAS_SIZE_X, CANVAS_SIZE_Y);
      canv[i]->SetLogx();
      canv[i]->SetLogy();
      hist_en[i]->Draw("colz");
      hist_en[i]->Write();
      
      TLatex label;
      label.SetNDC();
      label.SetTextSize(0.05);
      label.SetTextAlign(30);
      label.DrawLatex(0.92,0.875, entry[i].c_str());
      label.SetTextAlign(11);
      
      float eff = hist_en[i]->Integral(hist_en[i]->FindBin(25),
				       hist_en[i]->GetNbinsX()) / hist_en[i]->GetEntries();
      
      //float eff_err = TMath::Sqrt(eff*(1-eff)/hist_en[i]->GetEntries());
      eff*=100;
      //eff_err*=100;
      label.DrawLatex(0.20,0.225, Form("#splitline{#epsilon=%4.1f%%}"
				       "{Mean=%3.1f#pm%3.1ffC}",
				       eff,
				       hist_en[i]->GetMean(),
				       hist_en[i]->GetMeanError()));

      //canv[i]->Modified();
      canv[i]->Print(Form("Energy_Plots/energy_PS_%s.png", channels[i].name.c_str()));
      canv[i]->Print(Form("Energy_Plots/energy_PS_%s.pdf", channels[i].name.c_str()));
      canv[i]->Print(Form("Energy_Plots/energy_PS_%s.C", channels[i].name.c_str()));
      delete canv[i];
    }
    
    // ********************************************************************************
    // ******** STANDARD (EVEN SPACED) BINNING ENERGY HISTOGRAM *******
    // ********************************************************************************
    if (bins) {
      canv_bins[i] = new TCanvas(Form("en_bin_canv_%s", channels[i].name.c_str()), "", CANVAS_SIZE_X, CANVAS_SIZE_X);
      canv_bins[i]->SetLogx(0);
      canv_bins[i]->SetLogy();
      hist_en_bins[i]->Draw("colz");
      hist_en_bins[i]->Write();
      
      TLatex label_bins;
      label_bins.SetNDC();
      label_bins.SetTextSize(0.05);
      label_bins.SetTextAlign(30);
      label_bins.DrawLatex(0.92, 0.875, Form( "%s bins", entry[i].c_str()));
      label_bins.SetTextAlign(11);
      
      float eff_bins = hist_en_bins[i]->Integral(hist_en_bins[i]->FindBin(25),
						 hist_en_bins[i]->GetNbinsX()) / hist_en_bins[i]->GetEntries();
      
      //float eff_err = TMath::Sqrt(eff*(1-eff)/hist_en[i]->GetEntries());
      eff_bins*=100;
      //eff_err*=100;
      label_bins.DrawLatex(0.20,0.225, Form("#splitline{#epsilon=%4.1f%%}"
					    "{Mean=%3.1f#pm%3.1ffC}",
					    eff_bins,
					    hist_en_bins[i]->GetMean(),
					    hist_en_bins[i]->GetMeanError()));

      //canv_bins[i]->Modified();
      canv_bins[i]->Print(Form("Energy_Plots/energy_PS_bins_%s.png", TString(channels[i].name.c_str()).ReplaceAll("-","_").Data()));
      canv_bins[i]->Print(Form("Energy_Plots/energy_PS_bins_%s.pdf", TString(channels[i].name.c_str()).ReplaceAll("-","_").Data()));
      canv_bins[i]->Print(Form("Energy_Plots/energy_PS_bins_%s.C", TString(channels[i].name.c_str()).ReplaceAll("-","_").Data()));
      delete canv_bins[i];
    }

    // ********************************************************************************
    // ******** PEDESTAL HISTOGRAM *******
    // ********************************************************************************
    if (pedestal) {
      canv_ped[i] = new TCanvas(TString(channels[i].name.c_str()).ReplaceAll("-","_").Data(), "", CANVAS_SIZE_X, CANVAS_SIZE_Y);
      canv_ped[i]->SetLogx();
      canv_ped[i]->SetLogy();
      hist_en_ped[i]->Draw("colz");
      hist_en_ped[i]->Write();
      
      TLatex label_ped;
      label_ped.SetNDC();
      label_ped.SetTextSize(0.05);
      label_ped.SetTextAlign(30);
      label_ped.DrawLatex(0.92, 0.875, Form( "%s ped", entry[i].c_str()));
      label_ped.SetTextAlign(11);
      
      // note the +1 in GetNbinsX to include overflow bin
      float fakerate_ped = hist_en_ped[i]->Integral(hist_en_ped[i]->FindBin(25),
						    hist_en_ped[i]->GetNbinsX()+1) / hist_en[i]->GetEntries();
      
      //float fakerate_err = TMath::Sqrt(eff*(1-eff)/hist_en[i]->GetEntries());
      fakerate_ped*=100;
      //fakerate_err*=100;
      label_ped.DrawLatex(0.20,0.225, Form("#splitline{fake rate=%4.1f%%}"
					   "{Mean=%3.1f#pm%3.1ffC}",
					   fakerate_ped,
					   hist_en_ped[i]->GetMean(),
					   hist_en_ped[i]->GetMeanError()));
      
      //canv_ped[i]->Modified();
      canv_ped[i]->Print(Form("Energy_Plots/energy_PS_ped_%s.png", TString(channels[i].name.c_str()).ReplaceAll("-","_").Data()));
      canv_ped[i]->Print(Form("Energy_Plots/energy_PS_ped_%s.pdf", TString(channels[i].name.c_str()).ReplaceAll("-","_").Data()));
      canv_ped[i]->Print(Form("Energy_Plots/energy_PS_ped_%s.C", TString(channels[i].name.c_str()).ReplaceAll("-","_").Data()));
      delete canv_ped[i];
    } // if (ped)

  } // channels loop

  // ********************************************************************************
  // ******** OVERLAY THE PEDESTALS *******
  // ********************************************************************************
  if (overlay) {

    TLegend* leg_allped = new TLegend(0.6,0.56,0.87,0.9,"","brNDC");
    leg_allped->SetTextSize(.04);

    TCanvas* canv_allped = new TCanvas("overlayed_pedestal", "", CANVAS_SIZE_X, CANVAS_SIZE_Y);
    canv_allped->SetLogx();
    
    float max_value = -99;
    int max_index = 0;
    for (unsigned int k = 0; k < NUMCHAN; k++) {
      hist_en_ped[k]->SetLineColor(color[k]);
      hist_en_ped[k]->Draw("same");
      leg_allped->AddEntry(hist_en_ped[k], entry[k].c_str(), "l");
      if (hist_en_ped[k]->GetMaximum()>max_value) {
	max_index = k;
	max_value = hist_en_ped[k]->GetMaximum();
      }
    }
    
    hist_en_ped[max_index]->Draw();
    for (unsigned int k = 0; k < NUMCHAN; k++)
      hist_en_ped[k]->Draw("same");
    leg_allped->Draw("same");    
    hist_en_ped[max_index]->Draw("same,axis");

    canv_allped->Print("Overlayed_Plots/energyPS_all_ped.png");
    canv_allped->Print("Overlayed_Plots/energyPS_all_ped.pdf");
    canv_allped->Print("Overlayed_Plots/energyPS_all_ped.C");
    delete canv_allped;
  }

  energy_hists->Close();

}

////////////////////////////////////////////////////////////////////////////////
// Time-slice plot
////////////////////////////////////////////////////////////////////////////////
void doTimeSlice(bool debug, const char* dir) {

  // Let us define the histogram and book them
  TH1F *hist_ts[NUMCHAN];
  TH1F *hist_tsF[NUMCHAN];

  if (channels.size()!= NUMCHAN) {
    cout << "Argh, something wrong!" << endl;
    return;
  }

  for (unsigned int i = 0; i < channels.size(); ++i) {
    hist_ts[i] = new TH1F( Form("ts_%s", channels[i].name.c_str()), "", 10, 0.5, 10.5);
    hist_tsF[i] = new TH1F( Form("tsF_%s", channels[i].name.c_str()), "", 10, 0.5, 10.5);
  }

  // Now we get the data
  TChain* chain = new TChain("slim");
  chain->Add(Form("%s/*_slim.root", dir));

  // Get the branches I need:
  vector<double> *xa = 0, *xc = 0, *ya = 0, *yc = 0;
  chain->SetBranchAddress("xA",&xa);
  chain->SetBranchAddress("xC",&xc);
  chain->SetBranchAddress("yA",&ya);
  chain->SetBranchAddress("yC",&yc);

  // Alberto: this was double in James's macro
  float pulse[NCH][NTS], ped[NCH];
  chain->SetBranchAddress("pulse",&pulse);
  chain->SetBranchAddress("ped",&ped);

  double intercept_X, slope_X;
  double intercept_Y, slope_Y;
  chain->SetBranchAddress("interceptX",&intercept_X);
  chain->SetBranchAddress("slopeX",&slope_X);
  chain->SetBranchAddress("interceptY",&intercept_Y);
  chain->SetBranchAddress("slopeY",&slope_Y);

  int run;
  chain->SetBranchAddress("run", &run);

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
    for (unsigned int i = 0; i < channels.size(); ++i) {
      // This ensures that the finger tiles are using the correct runs
      if (run < 3410 && i >= 3 && i <= 6)
        continue;

      // Use the time slices numbers contained in TIMESLICES
      // and remove pedestal times # of time slices used 
      double energy_ps = 0;
      for (auto ts : TIMESLICES)
        energy_ps+=pulse[channels[i].chan][ts];

      // Here we remove the pedestal: there is a single number, corresponding to the
      // pedestal in one time slice; hence, we multiply it by the number of time slices
      // that form our signal pulse (defined in TIMESLICES)
      energy_ps-=TIMESLICES.size()*ped[channels[i].chan];

      // Very important: we require the hit to be fiducial
      if (!isFiducial( i, x_hit, y_hit))
        continue;

      // The fill here puts the pulse - ped in a specific bin
      for (int t=0;t<NTS;++t) {
	// i is index on channels: <3 and 7 are sigma tiles; 3-6 are finger tiles
        if (i < 3 || i == 7) {
          hist_ts[i]->Fill(t+1,
                           pulse[channels[i].chan][t]-
                           ped[channels[i].chan]);
        }
        else {
          hist_tsF[i]->Fill(t+1,
                            pulse[channels[i].chan][t]-
                            ped[channels[i].chan]);
        }
      }
    } // loop on channels
  } // loop on events

  // Let us also get the histogram with max value...
  int max_index = 0;
  int max_indexF = 0;
  float max_value = -99;
  float max_valueF = -99;
    
  TLegend* leg = new TLegend(0.2,0.7,0.4,0.9,"","brNDC");
  TLegend* leg_fing = new TLegend(0.2,0.7,0.4,0.9,"","brNDC");
  for (unsigned int i = 0; i < channels.size(); ++i) {
    if ((i < 3) || (i == 7)) {
      hist_ts[i]->SetLineWidth(2);
      hist_ts[i]->SetLineColor(color[i]);
      hist_ts[i]->SetLineStyle(style[i]);
	
      hist_ts[i]->GetXaxis()->SetTitle("Time Slice [25ns]");
      hist_ts[i]->GetXaxis()->SetNdivisions(NTS);
      hist_ts[i]->GetYaxis()->SetTitle("Charge [fC]");
      hist_ts[i]->GetYaxis()->SetTitleOffset(1.4);
      hist_ts[i]->GetYaxis()->SetMaxDigits(5);
      leg->AddEntry(hist_ts[i],entry[i].c_str(),"l");
	
      if (hist_ts[i]->GetMaximum()>max_value) {
	max_index = i;
	max_value = hist_ts[i]->GetMaximum();
      }
    }
    else { // For finger tiles
      hist_tsF[i]->SetLineWidth(2);
      hist_tsF[i]->SetLineColor(color[i]);
      hist_tsF[i]->SetLineStyle(style[i]);
	
      hist_tsF[i]->GetXaxis()->SetTitle("Time Slice [25ns]");
      hist_tsF[i]->GetXaxis()->SetNdivisions(NTS);
      hist_tsF[i]->GetYaxis()->SetTitle("Charge [fC]");
      hist_tsF[i]->GetYaxis()->SetTitleOffset(1.4);
      hist_tsF[i]->GetYaxis()->SetMaxDigits(4);
      leg_fing->AddEntry(hist_tsF[i],entry[i].c_str(),"l");
	
      if (hist_tsF[i]->GetMaximum()>max_valueF) {
	max_indexF = i;
	max_valueF = hist_tsF[i]->GetMaximum();
      }
    }
  }

  TCanvas* canv_ts = new TCanvas( "canv_ts", "", CANVAS_SIZE_X, CANVAS_SIZE_Y);    
  leg->SetTextSize(0.05);
  hist_ts[max_index]->Draw("hist");
  for (unsigned int i=0;i<channels.size(); i++) {
    if ((i < 3) || (i == 7))
      hist_ts[i]->Draw("hist,same");
  }
  leg->Draw("same");
  hist_ts[max_index]->Draw("axis,same");
    
  canv_ts->Print("Time_Slice_Plots/ts.png");
  canv_ts->Print("Time_Slice_Plots/ts.pdf");
  canv_ts->Print("Time_Slice_Plots/ts.C");
  delete canv_ts;    
    
  TCanvas* canv_tsF = new TCanvas("canv_tsF","",CANVAS_SIZE_X,CANVAS_SIZE_Y);
  leg_fing->SetTextSize(0.05);
  hist_tsF[max_indexF]->Draw("hist");
  for (unsigned int i=0;i<channels.size(); i++) {
    if (i > 2 || i < 7)
      hist_tsF[i]->Draw("hist,same");
  }
  leg_fing->Draw("same");
  hist_tsF[max_indexF]->Draw("axis,same");
    
  canv_tsF->Print("Time_Slice_Plots/tsF.png");
  canv_tsF->Print("Time_Slice_Plots/tsF.pdf");
  canv_tsF->Print("Time_Slice_Plots/tsF.C");
  delete canv_tsF;

  return;
}
