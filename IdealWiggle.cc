/*Plots an ideal wiggle plot, fits a function, calculates residual, performs FFT of the residual , derived from Talal/Becky's wiggle.cc All units are in us!
  SG 2018 */

#include <iostream>
#include "TH1F.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TVirtualFFT.h"

using namespace std;


//Wiggle fit function
Double_t fit(Double_t *v, Double_t *par) {
  return par[0]*TMath::Exp(v[0]*par[1]) * (1+par[2]*TMath::Cos(par[3]*v[0]+par[4]));
 }
//Exponential fit function
Double_t fitExp(Double_t *v, Double_t *par) {
  return par[0]*TMath::Exp(v[0]*par[1]);
 }

// Wiggle function - can inject a small signal 
Double_t Wiggle(Double_t *x, Double_t *par){
  
  Double_t time     = x[0];             // Leave time values free
  double tau       = 2.2;              // Lifetime of the muon at rest (us)
  double gamma     = 29.3;              // Magic gamma
  double A         = 0.05;              // Amplitude
  double omega     = 2 * TMath::Pi() / 4.2;   // Time of a single wiggle ~ 4.2 us
  double omega2     = 2 * TMath::Pi() / 0.2; // Random additional frequency
  double phase     = TMath::Pi()/2;     //Phase angle 
  double N = par[0];                        

  //    Double_t Npositrons = N * exp( - time / (tau * gamma)) * (1 + A * cos ( (omega * time) + phase)) ;
  Double_t Npositrons =  N * exp(- time / (tau * gamma) ) * (1 + A * (cos ( (omega * time) + phase ) +  0.1 * (cos ( (omega2 * time) + phase))));
 return Npositrons;

}

int main() {
  // 30x wiggles (900 spin precessions)
  const double totalTime = 30.0 * 4.2;

  TF1 *f1 = new TF1("f1", Wiggle , 0 , totalTime , 1); // Book ideal wiggle plot
  f1 -> SetNpx(10000); // Prevent undersampling
  f1 -> SetParameter(0,1); // You have to force the normalisation to be one.

  int nBins = 1000;
  
  TH1D *h1 = new TH1D("h1", "hist" , nBins , 0 , totalTime); // Book ideal histogram
  TH1D *h2 = new TH1D("h2", "hist" , nBins , 0 , totalTime); // Book residual

  //Instigate pseudorandom number generator 
  int seed = 12345; 
  TRandom3 *rand = new TRandom3(seed);
 
  // Loop over "events"
  double numEvents = 1e+8;

  for(int i=0; i < numEvents; i++) {
    double r1 = f1 -> GetRandom(0, totalTime);
    h1 -> Fill(r1);
  }  
  
  double functionIntegral = f1 -> Integral(0, totalTime);
  h1 -> Scale(1./functionIntegral);
  cout << "function integral " << functionIntegral << endl;
  
  //Perform the fits 
  TF1 *fit1 = new TF1("fit1", fit, 0, totalTime, 5);
  // Usually uneeded exponetial fit
  TF1 *fit2 = new TF1("fit2", fitExp, 0, totalTime, 2);

  fit1 -> SetNpx(10000);
  fit1 -> SetNpx(10000);
  fit1 -> SetLineWidth(1);
  fit2 -> SetLineWidth(1);
  fit2 -> SetLineColor(kGreen+2);

  //Fit parameters
  if (numEvents == 1e+7) {
   
    fit1 -> SetParLimits(3,1.4,1.5);
  }
  else if (numEvents == 1e+8) {
    // If pure
    fit1 -> SetParLimits(0, 3.5e3, 4.5e3);
    //Fit parameters 1e+8 events if 1 MHz injected
    // fit1 -> SetParLimits(2,0.040,0.050);
   fit1 -> SetParLimits(3,1.4,1.5);
  }
  else { 
    cout << "Fit parameters for this number of events do not exists, sorry.  SG"<<endl;
  }

  h1 -> Fit(fit1);
  
  //Get the residual                                          
  for (int i = 1; i < nBins; i++){

      double time = h1 -> GetBinCenter(i);
      double measured = h1 -> GetBinContent(i);
      double error = h1 -> GetBinError(i);
      double fitted = fit1 -> Eval(time);

      if (error > 0){
        h2 -> SetBinContent(i, measured-fitted); //residual SG   
      }
  }

  //get FFT of residual                                                       
  TH1 *hm = 0;
  TVirtualFFT::SetTransform(0);
  hm = h2 -> FFT(hm, "MAG");

  //Rescale x-axis by dividing by the function domain
              
  TAxis *xaxis = hm -> GetXaxis();
  double *ba = new double[nBins+1];
  xaxis -> GetLowEdge(ba);
  double Scale = 1./(totalTime - 0);
  ba[nBins] = ba[nBins-1] + xaxis -> GetBinWidth(nBins);

  for (int i = 0; i < nBins + 1; i++) {
       ba[i] *= Scale;
  }
 
  TH1F* fftResidual = new TH1F(hm -> GetName(), hm -> GetTitle(), nBins, ba);
  for (int i = 0; i <= nBins; i++) {
      fftResidual -> SetBinContent(i, hm -> GetBinContent(i));
       fftResidual -> SetBinError(i, hm -> GetBinError(i));
  }
  
  fftResidual -> SetTitle(";Frequency [MHz];Normalised Magnitude [Arb Units]");
  fftResidual -> SetStats(0);
  fftResidual -> SetName("residualFFT");
  fftResidual -> Scale(1.0 / fftResidual -> Integral());
 

  //Calculate Nyquist frequency, which is twice the highest frequeny in the signal or half of the sampling rate.                                                                                            
  //...the maximum frequency before sampling errors start              

  double binWidth = totalTime / nBins ;
  double sampleRate = 1 / binWidth;
  double nyquistFreq = 0.5 * sampleRate;
  cout << "binWidth " <<binWidth<<" us"<<endl;
  cout << "sampleRate " <<sampleRate<<" MHz"<<endl;
  cout << "nyquistFreq " <<nyquistFreq<<" MHz"<<endl;

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
  //fit2 -> Draw("SAME");
  //  c2 -> SetLogy();
  //  c2 -> SaveAs("figures/injected_1MHz_1e7/injectedFittedHistogram.eps");
   c2 -> SaveAs("figures/injected_2_1e8/injectedFittedHistogram.eps");
  //  c2 -> SaveAs("figures/pure_1e8/fittedHistogram.eps");
  
  TCanvas *c3 = new TCanvas();
  // gPad->SetGrid();
  h2 -> GetYaxis() -> SetMaxDigits(2);
  h2 -> SetTitle(";Time [#mus];Residual");
  h2 -> GetYaxis() -> SetMaxDigits(2);
  h2 -> SetStats(kFALSE);
  h2 -> Draw();                                                         
  //   c3 -> SaveAs("figures/injected_1MHz_1e7/residualInjected.eps");
   c3 -> SaveAs("figures/injected_2_1e8/residualInjected.eps");
  // c3 -> SaveAs("figures/pure_1e8/residual.eps");
  
  TCanvas *c4 = new TCanvas(); 
  fftResidual->SetTitle("FFT of Residual;Frequency [MHz];Magnitude");//, Spin Precession Frequency: #mbox{0.238} MHz;Frequency [MHz];Magnitude"); 
  //  fftResidual -> SetAxisRange(0,nyquistFreq,"X");		
  //c4 -> SetLogy();
  fftResidual -> GetYaxis() -> SetMaxDigits(2);
  fftResidual->GetXaxis()->SetMaxDigits(2);
  fftResidual -> SetMarkerColor(kRed+2);
  fftResidual -> SetLineColor(kRed+2);
  fftResidual -> Draw("HIST L SAME");
    c4 -> SaveAs("figures/injected_2_1e8/residualInjectedFFT.pdf");
  //  c4 -> SaveAs("figures/injected_1MHz_1e7/residualInjectedFFT.eps");
  // c4 -> SaveAs("figures/pure_1e8/residuaFFT.eps");

  delete f1;
  delete h1;
  delete h2;
  delete c1;
  delete c2;
  delete c3;
  delete hm;
  delete fftResidual;
  delete fit1;
  delete fit2;
  return 0;  
}
