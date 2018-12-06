//A place to test things, derived from Talal/Becky's wiggle.cc 
//All units are in ns
//SG 2018

#include <iostream>
#include "TH1F.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TVirtualFFT.h"

using namespace std;

Double_t Wiggle(Double_t *x, Double_t *par){
  
  Double_t time     = x[0];             // Leave time values free
  double tau       = 2200;              // Lifetime of the muon at rest (ns)
  double gamma     = 29.3;              // Magic gamma
  double A         = 0.05;              // Amplitude
  double omega     = 2 * TMath::Pi() / 4200;   // Time of a single wiggle ~ 4200 ns (omega=2pi/t)
  double omega2     = 2 * TMath::Pi() / 5000; // Random additional high frequency with T = 1000 ns
  double omega3     = 2 * TMath::Pi() / 500; // Random additional high frequency with T = 10000 ns
  double phase     = TMath::Pi()/2;     //Phase angle 
  double N = par[0];                        

   Double_t Npositrons = N * exp( - time / (tau * gamma)) * (1 + A * cos ( (omega * time) + phase)) ;
  // Double_t Npositrons =  N * exp(- time / (tau * gamma) ) * (1 + A * (cos ( (omega * time) + phase ) + cos ( (omega2 * time) + phase) + cos( (omega3 * time) + phase)));
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
  double numEvents = 1e+8;
  for(int i=0; i<numEvents;i++) {
    double r1 = f1 -> GetRandom(0,totalTime);
    h1 -> Fill(r1);
  }  
  
  double functionIntegral = f1->Integral(0,totalTime);
  h1->Scale(1./functionIntegral);
  cout<<functionIntegral<<endl;

  //get FFT of residual
  TH1 *hm =0;
  TVirtualFFT::SetTransform(0);
  hm=h1->FFT(hm,"MAG"); 

  //Rescale x-axis by dividing by the function domain 
  TAxis* xaxis = hm->GetXaxis();
  double* ba = new double[nBins+1];
  xaxis->GetLowEdge(ba);
  double Scale = 1./(totalTime - 0);
  ba[nBins] = ba[nBins-1] + xaxis->GetBinWidth(nBins);
  
  for (int i = 0; i < nBins+1; i++) {
       ba[i] *= Scale;
  }

  TH1F* fftResidual = new TH1F(hm->GetName(), hm->GetTitle(), nBins, ba);
  for (int i = 0; i <= nBins; i++) {
      fftResidual->SetBinContent(i, hm->GetBinContent(i));
       fftResidual->SetBinError(i, hm->GetBinError(i));
  }
 
  fftResidual->SetTitle(";Frequency [Hz];Normalised Magnitude [Arb Units]");
  fftResidual->SetStats(0);
  fftResidual->SetName("residualFFT");
  fftResidual->Scale(1.0 / fftResidual->Integral());

  TCanvas *c1 = new TCanvas();
  gPad->SetGrid();
  f1 -> SetTitle("Ideal Wiggle Plot;Time[ns];Log Normalised Units");
  f1 -> Draw();
  //  c1 -> SetLogy();
  c1 -> SaveAs("figures/IdealWiggle.pdf");


  TCanvas *c2 = new TCanvas();
  // gPad->SetGrid();
  h1 -> SetTitle(";Time [ns];Normalised Units");
  h1 -> GetYaxis() -> SetMaxDigits(2);
  h1 -> Draw();
  //  c2 -> SetLogy();
  c2 -> SaveAs("figures/IdealHistogram.pdf");

  TCanvas *c3 = new TCanvas(); 
  //gPad->SetGrid();
  //fftResidual -> SetTitle("FFT, Injected: #mbox{2#times10^{-3}} GHz & #mbox{2#times10^{-4}} GHz, g-2 freq: #mbox{2.38#times10^{-4}} GHz;Frequency [GHz];Magnitude [norm. units]");
  fftResidual->SetTitle("FFT, g-2 freq: #mbox{2.38#times10^{-4}} GHz;Frequency [GHz];Magitude [norm. units]"); 
  fftResidual -> SetAxisRange(0,0.004,"X");		
  c3 -> SetLogy();
  fftResidual->GetXaxis()->SetMaxDigits(2);
  fftResidual -> SetMarkerColor(kRed+2);
  fftResidual -> SetLineColor(kRed+2);
  fftResidual -> Draw("HIST L SAME");
  c3 -> SaveAs("figures/FFT.pdf");

  delete f1;
  delete h1;
  delete c1;
  delete c2;
  delete c3;
  delete hm;
  delete fftResidual;
  return 0;  
}
