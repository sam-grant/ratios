//Fourier transform beat pattern
#include <iostream>
#include "TH1F.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TVirtualFFT.h"
#include "TRandom3.h"

using namespace std;
Double_t Total_Time = 30.0 * 4200.0;

Double_t CosBasic(Double_t *x, Double_t *par){

  Double_t time     = x[0];
  double tau       = 2200;              // Lifetime of the muon at rest (ns)
  double gamma     = 29.3;              // Magic gamma
  double A         = 0.05;              // Amplitude
  double omega     = 2 * TMath::Pi() / 4200;   // Time of a single wiggle ~ 4200 ns (omega=2pi/t)
  double phase     = TMath::Pi()/2;
  double tau_gamma = tau * gamma;       //Making life easier
  double N = 1; 

  return  cos (time);
  
}

int main() {
   //prepare the canvas for drawing
  int seed = 12345; //SPG
   TRandom3 *rand = new TRandom3(seed);
   //add function
TF1 *f1 = new TF1("f1", CosBasic, 0, Total_Time, 1); 

    
    TH1D *h1 = new TH1D("h1", "Fourier Transform", 10000, 0, Total_Time);
    //  double x;
 
    //Fill the histogram with function values
    for (Int_t i=0; i<100000; i++){
      double x = f1->GetRandom();
	h1->Fill(x);//hsin->SetBinContent(i+1, f1->Eval(CosBasic));
  }

   TH1 *hm;
   //SetTransform(0);
   hm = h1->FFT(hm, "MAG");

    //TVirtualFFT::SetTransform(0);

  TCanvas* c1 = new TCanvas();
  gPad->SetGrid();
  h1 -> SetTitle("Hist cos (time)");
  h1 -> GetXaxis() -> SetTitle("Time (ns)");
  h1 -> GetYaxis() -> SetTitle("Number of Positrons N");
  //f1 -> GetYaxis() -> SetRangeUser(0, 2E11* 1.1);
  h1 -> Draw("C");
  c1 -> SaveAs("figures/CosBasicHist.pdf");

 TCanvas* c2 = new TCanvas();
  gPad->SetGrid();
  // hm->SetTitle("Magnitude of the 1st transform");
  // hm->Draw()
  hm -> SetTitle("FFT cos (time)");
  hm -> GetXaxis() -> SetTitle("Time (ns)");
  hm -> GetYaxis() -> SetTitle("Number of Positrons N");
  //f1 -> GetYaxis() -> SetRangeUser(0, 2E11* 1.1);
  hm -> Draw();
  c2 -> SaveAs("figures/CosBasicFFT.pdf");

  // delete f7;
  return 0;
}
