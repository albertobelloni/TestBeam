/*//////////////////////////////////////////////////////////////////////////////

Simple fit function for SiPM spectra
Uses TSpectrum to find peaks, and fits a Gaussian on top of each peak

Example usage:

 [] .L peakfinder.C+
 [] peakfinder("energy_hists.root","energy_tree_SCSN_81F1")

One can also run all the fits with one command:

 [] .L peakfinder.C+
 [] run_peakfinder("energy_hists.root")

INITIAL IDEA:

The shape parameters of each Gaussian are related, but the normalizations
are independent:

  mean_i = mean_1 + i * gain
  sigma_i² = sigma_0² + i * sigma_1²

The pedestal is set to be a Gaussian centered in mean_ped, with width sigma_ped

Parameters:

  0:    pedestal normalization
  1:    pedestal mean
  2:    pedestal sigma
  3:    gain
  4:    sigma_0
  5:    sigma_1
  6-15: normalization of first 10 peaks

WHAT WAS DONE:

The fit function is a plain sum of Gaussians, the parameters are:

  3*p:   normalization of p-th peak
  3*p+1: mean of p-th peak
  3*p+2: width of p-th peak

 *//////////////////////////////////////////////////////////////////////////////

#include "TCanvas.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TLine.h"
#include "TLatex.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include "TFitResult.h"
#include "TFile.h"
#include "TMath.h"
#include "string.h"
#include <stdio.h>

#define DEBUG false

using std::cout;
using std::endl;

int NPEAKS = 20;
Double_t fpeaks(Double_t *x, Double_t *par) {
  Double_t result = 0.;
  for (Int_t p=0;p<NPEAKS;p++) {
    Double_t norm  = par[3*p+0];
    Double_t mean  = par[3*p+1];
    Double_t sigma = par[3*p+2];
    result += norm*TMath::Gaus(x[0],mean,sigma,true);
  }
  return result;
}

// Helper function used to print the fit results
// for the DN-18-007 note and the HCAL DPG presentation
int idx(const char* treename);

// The main function: this is the one you call!
void peakfinder(const char* filename,
		const char* treename) {

  
  TFile *file = TFile::Open(filename);
  if (!file) {
    std::cout << "File " << filename << " does not seem to exist\n";
    return;
  }

  TTree *tree = (TTree*)file->Get(treename);
  if(!tree) {
    std::cout << "The tree pointer is null: is the tree name correct?\n";
    return;
  }
  if(strchr(tree->ClassName(),'H')) {
    std::cout << "The tree is a TH1F: are both the fit type and "
	      << "the tree name correct?\n";
      return;
  }
  
  TCanvas *canvas = new TCanvas("canvas","",800,600);
  canvas->Draw();

  const char* histname = Form("hist_%s",treename);
  TH1F *hist = new TH1F(histname,"",100,-50,550);
  tree->Project(histname,"x"); // x is the variable inside the relevant tree

  TSpectrum *spec = new TSpectrum(2*NPEAKS);

  // NOTE: using sigma=1 to search for peaks, but then use sigma=10
  // as a fit function parameter: it is much closer to the actual
  // fit result (sigma between 7 and 9)
  int found = spec->Search(hist,1,"nobackground",0.001);
  printf("\n\nFound %d peaks\n\n",found);

  // Let us take the initial parameters for the fit from found peaks
  double parms[3*found];
  NPEAKS=0; // why? why?!?
  
  double *xpeaks = spec->GetPositionX();
  for (int p=0;p<found;p++) {
    double xp = xpeaks[p];
    double yp = hist->GetBinContent(hist->GetXaxis()->FindBin(xp));
    if (xp<-10) continue;
    parms[3*NPEAKS+0] = 0.8*yp*sqrt(2*TMath::Pi())*8;
    parms[3*NPEAKS+1] = xp;
    parms[3*NPEAKS+2] = 8;
    NPEAKS++;
  }
  printf("Found %d useful peaks to fit\n\n",NPEAKS);

  if(DEBUG) {
    printf("Peak\tNormalization\tMean\tWidth\n");
    for (int p=0;p<NPEAKS;p++)
      printf("%d\t%12.1f\t%.1f\t%.1f\n",p,parms[3*p],parms[3*p+1],parms[3*p+2]);
    printf("\n");
  }

  // Now we build the fit function
  double fit_range_max = 375;
  if (strcmp(treename,"energy_tree_EJ_200")   ==0 ||
      strcmp(treename,"energy_tree_EJ_260")   ==0 ||
      strcmp(treename,"energy_tree_EJ_260_2P")==0 ||
      strcmp(treename,"energy_tree_SCSN_81S") ==0   )
    fit_range_max=500;

  if (strcmp(treename,"energy_tree_SCSN_81F3")==0)
    fit_range_max=400;
  if (strcmp(treename,"energy_tree_SCSN_81F4")==0)
    fit_range_max=395;

  TF1 *func = new TF1("func",fpeaks,-40,fit_range_max,3*NPEAKS);

  // We may have more than the default 25 parameters
  TVirtualFitter::Fitter(hist,3*NPEAKS);

  // The lines below do not seem to work
  //int exreturn = 0;
  //TVirtualFitter::GetFitter()->ExecuteCommand("SET STR 2",NULL,exreturn);

  func->SetParameters(parms);
  // Force normalization to be positive
  for (int p=0;p<NPEAKS;p++) {
    func->SetParLimits(3*p+0,0,2.5e6);

    func->SetParName(3*p+0,Form("Norm_%d",p));
    func->SetParName(3*p+1,Form("Mean_%d",p));
    func->SetParName(3*p+2,Form("Width_%d",p));
  }

  // Let us also fix to 0 the normalization of the Gaussians
  // which are out of the fitting range (i.e., if mean-2*width<fit_range_max)
  // I hope this will make the fit more sensible
  for (int p=0;p<NPEAKS;p++)
    if (parms[3*p+1]-2*parms[3*p+2]>fit_range_max) {
      func->FixParameter(3*p, 0);
      func->FixParameter(3*p+1, 500);
      func->FixParameter(3*p+2,  10);
    }

  // The following line makes the function smoother, when drawn
  func->SetNpx(1000);

  // THIS IS THE FIT:
  TFitResultPtr fit = hist->Fit("func","RS");

  // Repeat fit until successful...
  while (fit->CovMatrixStatus()!=3 ||
	 fit->Status()!=0 ||
	 fit->IsValid()==0)
    fit = hist->Fit("func","RS");

  // The trick of the experienced physicist: since MIGRAD seems to
  // print on the screen that the covariance matrix is approximate, but
  // still returns a perfect covariance matrix status, let us repeat
  // the fit by hand for those tiles we know need a small push
  if (strcmp(treename,"energy_tree_EJ_200")==0||
      strcmp(treename,"energy_tree_SCSN_81F1")==0||
      strcmp(treename,"energy_tree_SCSN_81F2")==0||
      strcmp(treename,"energy_tree_SCSN_81F3")==0)
    fit = hist->Fit("func","RS");

  // At this point, all fits are ok, with the exception of:
  // SCSN_81F4: uncertainty 1.8%
  // EJ_260: unceratinty 3.0%

  func->SetLineWidth(4);
  func->Draw("same");

  // Darn, I really need to do some sorting...
  int sorted[NPEAKS];
  double peak_loc[NPEAKS];
  for (int p=0;p<NPEAKS;++p) {
    sorted[p]=0;
    peak_loc[p]=func->GetParameter(3*p+1);
  }
  TMath::Sort(NPEAKS,peak_loc,sorted,false);

  // Let us draw one Gaussian for each peak, to check
  // that fit was successful
  for (int p=0;p<NPEAKS;++p) {
    if (DEBUG) {
      printf("\n sorted? p=%d, sorted[p]=%d, "
	     "peak_loc[p]=%f, peak_loc[sorted[p]]=%f, "
	     "norm[p]=%f, norm[sorted[p]]=%f, "
	     "mean[p]=%f, mean[sorted[p]]=%f\n",
	     p, sorted[p], peak_loc[p], peak_loc[sorted[p]],
	     func->GetParameter(3*p),
	     func->GetParameter(3*sorted[p]),
	     func->GetParameter(3*p+1),
	     func->GetParameter(3*sorted[p]+1));
    }
    auto peak_gauss = new TF1(Form("peak_gauss_%d",p),"[area] * "
		      "ROOT::Math::normal_pdf(x, [sigma], [mean]) ",
		      -50, 550);
    peak_gauss->SetParameter("area",func->GetParameter(3*p));
    peak_gauss->SetParameter("mean",func->GetParameter(3*p+1));
    peak_gauss->SetParameter("sigma",func->GetParameter(3*p+2));
    peak_gauss->SetLineColor(p % 9 + 1);
    peak_gauss->SetLineStyle(p >9?2:1);
    peak_gauss->SetNpx(1000);
    peak_gauss->Draw("same");
  }

  // Calculate the average number of photo-electrons
  // Need to skip the pedestal peak
  double pe_numerator = 0;
  double pe_denominator = 0;
  double gain_estimate = 0;
  for (int p=0;p<NPEAKS;++p) {
    // Skip the pedestal
    if (func->GetParameter(3*sorted[p]+1)<20) continue;

    // Skip adding peak if normalization higher than twice the previous one
    // it means we have a problem with fit
    // Usually also the peak after it are bad, let us break...
    if (func->GetParameter(3*sorted[p])>
	func->GetParameter(3*sorted[p-1])*2)
	break;

    // Skip adding peak if it is outside the fitting range
    // (it should be a rather small correction, even if the fit were good)
    if (func->GetParameter(3*sorted[p]+1)
	-2*func->GetParameter(3*sorted[p]+2)>fit_range_max)
      continue;

    pe_numerator += p * func->GetParameter(3*sorted[p]);
    pe_denominator += func->GetParameter(3*sorted[p]);

    if (DEBUG)
      printf("\n p=%d, peak_loc[sorted[p]]=%f, "
	     "mean[sorted[p]]=%f, "
	     "norm[sorted[p]]=%f\n",
	     p, peak_loc[sorted[p]],
	     func->GetParameter(3*sorted[p]+1),
	     func->GetParameter(3*sorted[p]));

    // Let us draw a line on top of each Gaussian used to calculate
    // average number of photo-electrons
    // The line goes from the minimum of the histogram to 30% more
    // of each Gaussian peak (or the maximum of the histogram,
    // whichever smaller)
    auto line_gauss = new TLine(func->GetParameter(3*sorted[p]+1),
				hist->GetMinimum(),
			        func->GetParameter(3*sorted[p]+1),
				hist->GetMaximum()<
				func->GetParameter(3*sorted[p]+0)/
				(TMath::Sqrt(2*TMath::Pi())*
				 func->GetParameter(3*sorted[p]+2))*1.3?
				hist->GetMaximum():
				func->GetParameter(3*sorted[p]+0)/
				(TMath::Sqrt(2*TMath::Pi())*
				 func->GetParameter(3*sorted[p]+2))*1.3);
    line_gauss->SetLineColor(sorted[p] % 9 + 1);
    line_gauss->SetLineStyle(sorted[p] >9?2:1);
    line_gauss->SetLineWidth(2);
    line_gauss->Draw("same");

    // Note: p=0 should correspond, after ordering the list of peaks,
    // to the pedestal peak: we should have skipped it already
    // (look at "continue" above: want Gaussian mean above 20fC)
    if (p>1)
      // add to gain estimate the difference between each pair of peaks
      gain_estimate +=
	func->GetParameter(3*sorted[p]+1)-
	func->GetParameter(3*sorted[p-1]+1);
  }
  gain_estimate/=(NPEAKS-1); // there are NPEAKS-1 pairs of neighboring peaks
  // Upon noting that the peaks are ordered by height, and not by location,
  // the algorithm above is sometimes wrong. Since the physics-driven fit
  // returns a gain between 41.6 and 39.9 for all tiles, let me use 41.0 here,
  // and accept an up to 3% mistake in the yield estimate
  gain_estimate = 41.;
	
  // It looks like we will need the following numbers a few times
  // Let us calculate them once and for all
  float avg_pe_fit = pe_denominator>0?pe_numerator/pe_denominator:-99;
  // Some machinery to calculate mean: want to use only the range with
  // at least 1 p.e.
  hist->GetXaxis()->SetRangeUser(25,600);
  float avg_pe_integ = gain_estimate>0?hist->GetMean()/gain_estimate:-99;
  hist->GetXaxis()->SetRangeUser(-50,550);

  printf("\n\nTile: %s\n",treename);
  printf("Average number of p.e. from fit: %.3f\n",avg_pe_fit);
  printf("Average number of p.e. from TH1::Mean: %.3f\n",avg_pe_integ);
  printf("Estimated gain: %.3f\n\n",gain_estimate);

  // Save formatted outputs to a file
  FILE *output;
  output = fopen("yield_results.txt","a");

  // Formatted output for presentation
  int row=idx(treename);
  fprintf(output,"PRESENTATION table.cell(%d,4).text = \"%.3f\"\n",
	  row,avg_pe_fit);
  fprintf(output,"PRESENTATION table.cell(%d,5).text = \"%.3f\"\n",
	  row,avg_pe_integ);

  // Formatted output for DN-18-007
  fprintf(output,
	  "DN-18-007_NOTE %d & %.3f & %.3f\n",row,avg_pe_fit,avg_pe_integ);

  fclose(output);

  // Minimal beautification of the histogram
  hist->GetXaxis()->SetTitle("Charge [fC]");
  hist->GetYaxis()->SetTitle("Events");
  hist->GetXaxis()->SetTitleOffset(1.1);
  hist->GetYaxis()->SetTitleOffset(1);
  hist->GetYaxis()->SetMaxDigits(3);
  hist->SetLineWidth(3);
  hist->Draw("hist,same");

  TLatex label;
  label.SetNDC();
  label.SetTextSize(0.05);
  label.SetTextAlign(31); // right-aligned; center-aligned
  label.DrawLatex(0.91,0.875,Form("<p.e.> hist: %.3f",avg_pe_integ));
  label.DrawLatex(0.91,0.825,Form("<p.e.> fit: %.3f",avg_pe_fit));
  label.DrawLatex(0.91,0.775,TString(treename).
		  ReplaceAll("energy_tree_","").
		  ReplaceAll("_","-"));

  canvas->Print(TString(treename).
		ReplaceAll("energy_tree","canvas_gaussfit_lin").
		Append(".root"));
  // NOTE: problems with C file...
  canvas->Print(TString(treename).
		ReplaceAll("energy_tree","canvas_gaussfit_lin").
		Append(".C"));
  canvas->Print(TString(treename).
		ReplaceAll("energy_tree","canvas_gaussfit_lin").
		Append(".pdf"));
  canvas->Print(TString(treename).
		ReplaceAll("energy_tree","canvas_gaussfit_lin").
		Append(".png"));

  canvas->SetLogy(true);
  canvas->Update();

  canvas->Print(TString(treename).
		ReplaceAll("energy_tree","canvas_gaussfit_log").
		Append(".root"));
  // NOTE: problems with C file...
  canvas->Print(TString(treename).
		ReplaceAll("energy_tree","canvas_gaussfit_log").
		Append(".C"));
  canvas->Print(TString(treename).
		ReplaceAll("energy_tree","canvas_gaussfit_log").
		Append(".pdf"));
  canvas->Print(TString(treename).
		ReplaceAll("energy_tree","canvas_gaussfit_log").
		Append(".png"));


  return;

}


// When you need to run all fits one after the other...
void run_peakfinder(const char* filename) {

  const char* treenames[8] =
    {"energy_tree_EJ_200",
     "energy_tree_EJ_260",
     "energy_tree_EJ_260_2P",
     "energy_tree_SCSN_81F1",
     "energy_tree_SCSN_81F2",
     "energy_tree_SCSN_81F3",
     "energy_tree_SCSN_81F4",
     "energy_tree_SCSN_81S"};

  for (int i=0;i<8;peakfinder(filename,treenames[i++]));

  return;

}

int idx(const char* treename) {

  const char* treenames[8] =
    {"energy_tree_EJ_200",
     "energy_tree_EJ_260",
     "energy_tree_EJ_260_2P",
     "energy_tree_SCSN_81F1",
     "energy_tree_SCSN_81F2",
     "energy_tree_SCSN_81F3",
     "energy_tree_SCSN_81F4",
     "energy_tree_SCSN_81S"};

  int i = 0;
  for (i=0;i<8;++i)
    if (strcmp(treename,treenames[i])==0) break;

  return i+1;

}
