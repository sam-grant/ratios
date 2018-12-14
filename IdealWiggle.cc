//A place to test things, derived from Talal/Becky's wiggle.cc 
//All units are in us!
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

Double_t fit(Double_t *v, Double_t *par) {
  return par[0]*TMath::Exp(v[0]*par[1]) * (1+par[2]*TMath::Cos(par[3]*v[0]+par[4]));
 }
					     
							  
Double_t Wiggle(Double_t *x, Double_t *par){
  
  Double_t time     = x[0];             // Leave time values free
  double tau       = 2.2;              // Lifetime of the muon at rest (ns)
  double gamma     = 29.3;              // Magic gamma
  double A         = 0.05;              // Amplitude
  double omega     = 2 * TMath::Pi() / 4.2;   // Time of a single wiggle ~ 4200 ns (omega=2pi/t)
  double omega2     = 2 * TMath::Pi() / 5; // Random additional high frequency with T = 1000 ns
  double omega3     = 2 * TMath::Pi() / 0.5; // Random additional high frequency with T = 10000 ns
  double phase     = TMath::Pi()/2;     //Phase angle 
  double N = par[0];                        

  Double_t Npositrons = N * exp( - time / (tau * gamma)) * (1 + A * cos ( (omega * time) + phase)) ;
  //   Double_t Npositrons =  N * exp(- time / (tau * gamma) ) * (1 + A * (cos ( (omega * time) + phase ) + cos ( (omega2 * time) + phase) + cos( (omega3 * time) + phase)));
 return Npositrons;

}

int main() {

  //For h1, setting nBins to 1000 and numEvents to 10000 creates an ok histogram

  const double totalTime = 30.0 * 4.2;

  TF1 *f1 = new TF1("f1", Wiggle , 0 , totalTime , 1); // Book ideal wiggle plot
  f1 -> SetNpx(10000); // Prevent undersampling
  f1 -> SetParameter(0,1); // You have to force the normalisation to be one.

  int nBins = 1000;
  TH1D *h1 = new TH1D("h1", "hist" , nBins , 0 , totalTime); // Book ideal histogram

  //Instigate pseudorandom number generator 
  int seed = 12345; 
  TRandom3 *rand = new TRandom3(seed);
 
  // Loop over "events"
  double numEvents = 1e+7;
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
 
  fftResidual->SetTitle(";Frequency [MHz];Normalised Magnitude [Arb Units]");
  fftResidual->SetStats(0);
  fftResidual->SetName("residualFFT");
  fftResidual->Scale(1.0 / fftResidual->Integral());


  //Calculate Nyquist frequency, which is twice the highest frequeny in the signal or half of the sampling rate.
  //...the maximum frequency before sampling errors start 

  double binWidth = totalTime / nBins ;
  double sampleRate = 1 / binWidth;
  double nyquistFreq = 0.5 * sampleRate;
  cout << "binWidth " <<binWidth<<" us"<<endl;
  cout << "sampleRate " <<sampleRate<<" per us"<<endl;
  cout << "nyquistFreq " <<nyquistFreq<<" MHz"<<endl;

  
  //Perform the fit 
  TF1 *fit1 = new TF1("fit1",fit,0,totalTime,5);
  fit1 -> SetNpx(10000);
  fit1 -> SetLineWidth(1);
  // fit1 -> SetParLimits(0,330,340); // Expect 
  // fit1 -> SetParLimits(1,-1.55,-1.5); // Expect ~2.2
  //fit1 -> SetParLimits(2,0.2,0.3); // Expect ~
   fit1 -> SetParLimits(3,1.4,1.6);// Expect ~1.5 MHz
  //fit1 -> SetParLimits(4,1,2);// Expect ~1.5 MHz
  h1->Fit(fit1);
  // fit1 -> SetNpx(10000);

  
  //Get the residual                                          
  /* for (int ibin(1); ibin < BinNum; ibin++){
      double time = h1->GetBinCenter(ibin);
      double measured = h1->GetBinContent(ibin);
      double error = h1->GetBinError(ibin);
      double fitted = fit->Eval(time);

      //cout << "time: " << time << ", meas: " << measured << ", err: " << error <\
< ", fitted: " << fitted << "\n";                                                  

      // meas - fitted / error                                                     

      if (error > 0){
        h2 -> SetBinContent(ibin,(measured-fitted)); //residual SG   
  */
  TCanvas *c1 = new TCanvas();
  //  gPad->SetGrid();
  f1 -> SetTitle("Ideal Wiggle Plot;Time [#mus];Number of Positrons");
  f1 -> Draw();
  //  c1 -> SetLogy();
  c1 -> SaveAs("figures/IdealWiggle.eps");


  TCanvas *c2 = new TCanvas();
  // gPad->SetGrid();
  h1 -> SetTitle(";Time [#mus];Number of Positrons");
  h1 -> GetYaxis() -> SetMaxDigits(2);
  // h1 -> SetAxisRange(0,500,"Y");
  h1 -> SetStats(kFALSE);
  h1 -> Draw();
  fit1 -> Draw("SAME");
  //  c2 -> SetLogy();
  c2 -> SaveAs("figures/IdeaHistogram.eps");

  TCanvas *c3 = new TCanvas(); 
  //gPad->SetGrid();
  //fftResidual -> SetTitle("FFT, Injected: #mbox{2#times10^{-3}} GHz & #mbox{2#times10^{-4}} GHz, g-2 freq: #mbox{2.38#times10^{-4}} GHz;Frequency [GHz];Magnitude [norm. units]");
  fftResidual->SetTitle("FFT, Spin Precession Frequency: #mbox{0.238} MHz;Frequency [MHz];Magnitude"); 
  fftResidual -> SetAxisRange(0,nyquistFreq,"X");		
  c3 -> SetLogy();
  fftResidual->GetXaxis()->SetMaxDigits(2);
  fftResidual -> SetMarkerColor(kRed+2);
  fftResidual -> SetLineColor(kRed+2);
  fftResidual -> Draw("HIST L SAME");
  c3 -> SaveAs("figures/FFT.eps");

  delete f1;
  delete h1;
  delete c1;
  delete c2;
  delete c3;
  delete hm;
  delete fftResidual;
  return 0;  
}
