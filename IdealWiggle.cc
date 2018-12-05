//A place to test things, derived from Talal/Becky's wiggle.cc 
//All units are in ns
//SG 2018

#include <iostream>
#include "TH1F.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TVirtualFFT.h"
//#include "fftw3.h"
using namespace std;

//const Double_t Total_Time = 30.0 * 4200.0; // Global

Double_t Wiggle(Double_t *x, Double_t *par){
  
  Double_t time     = x[0];             // Leave time values free
  double tau       = 2200;              // Lifetime of the muon at rest (ns)
  double gamma     = 29.3;              // Magic gamma
  double A         = 0.05;              // Amplitude
  double omega     = 2 * TMath::Pi() / 4200;   // Time of a single wiggle ~ 4200 ns (omega=2pi/t)
  double omega2     = 2 * TMath::Pi() / 10; // Random additional high frequency with T = 1000 ns
  double omega3     = 2 * TMath::Pi() / 10000; // Random additional low frequency with T = 10000 ns
  double phase     = TMath::Pi()/2;     //Phase angle 
  double N = par[0];                        

  //  Double_t Npositrons = N * exp( - time / (tau * gamma)) * (1 + A * cos ( (omega * time) + phase)) ;
   Double_t Npositrons =  N * exp(- time / (tau * gamma) ) * (1 + A * (cos ( (omega * time) + phase ) + cos ( (omega2 * time) + phase) + cos( (omega3 * time) + phase)));
 return Npositrons;
}

int main() {

  //For h1, setting nBins to 1000 and numEvents to 10000 creates an ok histogram

  const double totalTime = 30.0 * 4200.0;

  TF1 *f1 = new TF1("f1", Wiggle , 0 , totalTime , 1); // Book ideal wiggle plot
  f1 -> SetNpx(10000); // Prevent undersampling
  f1 -> SetParameter(0,1); // You have to force the normalisation to be one.

  int nBins = 1000;
  TH1D *h1 = new TH1D("h1", "hist" , nBins , 0 , totalTime); // Book ideal histogram

  //Instigate pseudorandom number generator 
  int seed = 12345; 
  TRandom3 *rand = new TRandom3(seed);

  // Loop over "events"
  double numEvents = 2e+7;
  for(int i=0; i<numEvents; i++) {
    double r1 = f1 -> GetRandom(0,totalTime);
    // double log = TMath::Log(r1);
    h1 -> Fill(r1);
  }  
  
  double functionIntegral = f1->Integral(0,totalTime);
  cout<<functionIntegral<<endl;


  //get FFT of residual

  TH1 *fft = 0;
  TVirtualFFT::SetTransform(0);
  fft = h1 -> FFT(fft,"MAG");
  //  fft -> SetStats(kFALSE); // ??



  TCanvas *c1 = new TCanvas();
  gPad->SetGrid();
  f1 -> SetTitle("Ideal Wiggle Plot;Time[ns];Log Normalised Units");
  f1 -> Draw();
  //  c1 -> SetLogy();
  c1 -> SaveAs("figures/IdealWiggle.pdf");


  TCanvas *c2 = new TCanvas();
  gPad->SetGrid();
  h1 -> SetTitle("Ideal Histogram;Time [ns];N");
  h1 -> Draw();
  //  c2 -> SetLogy();
  c2 -> SaveAs("figures/IdealHistogram.pdf");

  TCanvas *c3 = new TCanvas(); 
  gPad->SetGrid();
  fft -> SetTitle("FFT;Frequency [Hz] ?;Magnitude ?");
  //  fft -> SetAxisRange(0,1/totalTime,"X");		
  //  fft -> SetAxisRange(0,1/sqrt(numEvents),"Y");
  c3 -> SetLogy();
  fft -> Draw();
  c3 -> SaveAs("figures/FFT.pdf");

  delete f1;
  delete h1;
  delete c1;
  delete c2;
  return 0;  
}
