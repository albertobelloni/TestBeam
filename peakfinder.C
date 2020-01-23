/*//////////////////////////////////////////////////////////////////////////////

Simple fit function for SiPM spectra
Uses TSpectrum to find peaks, and fits a Gaussian on top of each peak

Example usage:

 [] .L peakfinder.C+
 [] peakfinder("energy_hists.root","energy_tree_SCSN_81F1")

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
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include "TFile.h"

using std::cout;
using std::endl;

Int_t NPEAKS = 10;
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

void peakfinder(const char* filename,
		const char* tilename) {

  
  TFile *file = TFile::Open(filename);
  if (!file) {
    std::cout << "File " << filename << " does not seem to exist\n";
    return;
  }

  TTree *tree = (TTree*)file->Get(tilename);
  if(!tree) {
    std::cout << "The tree pointer is null: is the tree name correct?\n";
    return;
  }
  if(strchr(tree->ClassName(),'H')) {
    std::cout << "The tree is a TH1F: are both the fit type and "
	      << "the tree name correct?\n";
      return;
  }
  
  const char* histname = Form("hist_%s",tilename);
  TH1F *hist = new TH1F(histname,"",100,-50,550);
  tree->Project(histname,"x"); // x is the variable inside the relevant tree

  TSpectrum *spec = new TSpectrum(2*NPEAKS);

  // NOTE: using sigma=1 to search for peaks, but then use sigma=10
  // as a fit function parameter: it is much closer to the actual
  // fit result (sigma between 7 and 9)
  int found = spec->Search(hist,1,"nobackground",0.001);
  printf("\n\nFound %d peaks\n\n",found);

  // Let us take the initial parameters for the fit from found peaks
  double parms[3*NPEAKS];
  NPEAKS=0; // why? why?!?
  
  Double_t *xpeaks = spec->GetPositionX();
  for (int p=0;p<found;p++) {
    Double_t xp = xpeaks[p];
    Double_t yp = hist->GetBinContent(hist->GetXaxis()->FindBin(xp));
    if (xp<-10) continue;
    parms[3*NPEAKS+0] = yp*sqrt(2*TMath::Pi())*10.;
    parms[3*NPEAKS+1] = xp;
    parms[3*NPEAKS+2] = 10;
    NPEAKS++;
  }
  printf("Found %d useful peaks to fit\n\n",NPEAKS);

  printf("Peak\tNormalization\tMean\tWidth\n");
  for (int p=0;p<NPEAKS;p++)
    printf("%d\t%12.1f\t%.1f\t%.1f\n",p,parms[3*p],parms[3*p+1],parms[3*p+2]);
  printf("\n");

  // Now we build the fit function
  TF1 *func = new TF1("func",fpeaks,-40,400,3*NPEAKS);

  // We may have more than the default 25 parameters
  TVirtualFitter::Fitter(hist,3*NPEAKS);

  func->SetParameters(parms);
  // Force normalization and widths to be positive
  for (int p=0;p<NPEAKS;p++) {
    func->SetParLimits(3*p+0,0,2e5);
    //func->SetParLimits(3*p+2,0,1e2);

    func->SetParName(3*p+0,Form("Norm_%d",p));
    func->SetParName(3*p+1,Form("Mean_%d",p));
    func->SetParName(3*p+2,Form("Width_%d",p));
  }

  func->SetNpx(1000);
  hist->Fit("func","LER");

  // Calculate the average number of photo-electrons
  // Need to skip the first peak, corresponding to pedestal
  double pe_numerator = 0;
  double pe_denominator = 0;
  double gain_estimate = 0;
  for (int p=1;p<NPEAKS;p++) {
    pe_numerator += p * func->GetParameter(3*p);
    pe_denominator += func->GetParameter(3*p);
    if (p>1)
      // add to gain estimate the difference between each pair of peaks
      gain_estimate += func->GetParameter(3*p+1)-func->GetParameter(3*p-2);
  }
  gain_estimate/=(NPEAKS-1); // there are NPEAKS-1 pairs of neighboring peaks
	
  printf("\nAverage number of p.e. from fit: %f\n",pe_numerator/pe_denominator);

  // Some machinery to calculate mean: want to use only the range with
  // at least 1 p.e.
  hist->GetXaxis()->SetRangeUser(25,600);
  printf("\nAverage number of p.e. from TH1::Mean: %f\n",
	 hist->GetMean()/gain_estimate);
  printf("Estimated gain: %f\n\n",gain_estimate);
  hist->GetXaxis()->SetRangeUser(-50,550);

}
