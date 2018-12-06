/* Somewhere the test TVirtualFFT - SG 2018 */

#include <iostream>
#include <TVirtualFFT.h>
#include <TF1.h>
#include <TH1D.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TH1F.h>
//#include <TPad.h>
using namespace std;

int main() {
  /*
  TH1F* EuropaFitter::RescaleAxis(TH1* input, Double_t Scale) {
     int bins = input->GetNbinsX();
     TAxis* xaxis = input->GetXaxis();
     double* ba = new double[bins+1];
     xaxis->GetLowEdge(ba);
     ba[bins] = ba[bins-1] + xaxis->GetBinWidth(bins);
     for (int i = 0; i < bins+1; i++) {
       ba[i] *= Scale;
     }
     TH1F* out = new TH1F(input->GetName(), input->GetTitle(), bins, ba);
     for (int i = 0; i <= bins; i++) {
       out->SetBinContent(i, input->GetBinContent(i));
       out->SetBinError(i, input->GetBinError(i));
     }
     return out;
}
  */
//A function to sample
  double start = 0;
  double end = 4*TMath::Pi();
  TF1 *fsin = new TF1("fsin", "sin(x)", start, end);//4*TMath::Pi());
   fsin->Draw();
   
    Int_t n=100;
    TH1D *hsin = new TH1D("hsin", "hsin", n+1, start, end);//0, 4*TMath::Pi());
    Double_t x;
    
       //Fill the histogram with function values
       for (Int_t i=0; i<=n; i++){
            x = (Double_t(i)/n)*(4*TMath::Pi());
            hsin->SetBinContent(i+1, fsin->Eval(x));
         }

     TCanvas *c0 = new TCanvas;
      hsin->Draw();
      c0->SaveAs("figures/hsin.pdf");

    fsin->GetXaxis()->SetLabelSize(0.05);
    fsin->GetYaxis()->SetLabelSize(0.05);
   
    //Compute the transform and look at the magnitude of the output
    TH1 *hm =0;
    TVirtualFFT::SetTransform(0);
    //    hm = hsin->FFT(hm, "MAG");
    //Rescale x axis
    hm = hsin->FFT(hm,"MAG");

    // TH1F* fftResidualInit = SetupFFT(residual, fitStartTime_, fitEndTime_);
    //hm = fftResidualInit->FFT(hm,"MAG");
    //int bins = input->GetNbinsX();
     TAxis* xaxis = hm->GetXaxis();
     double* ba = new double[n+1];
     xaxis->GetLowEdge(ba);
     double Scale = 1./(end - start);
     ba[n] = ba[n-1] + xaxis->GetBinWidth(n);
     for (int i = 0; i < n+1; i++) {
       ba[i] *= Scale;
     }
     TH1F* fftResidual = new TH1F(hm->GetName(), hm->GetTitle(), n, ba);
     for (int i = 0; i <= n; i++) {
      fftResidual->SetBinContent(i, hm->GetBinContent(i));
      fftResidual->SetBinError(i, hm->GetBinError(i));
      }
     // TH1F* fftResidual = RescaleAxis(hm,1./(end - start));
    //TH1F* fftResidual = RescaleAxis(hm,1./(end - start));
    fftResidual->SetTitle(";Frequency (MHz);Magnitude [Arb Units]");
    fftResidual->SetStats(0);
    fftResidual->SetName("residualFFT");
    fftResidual->Scale(1.0 / fftResidual->Integral());
    // hm->SetTitle(";Frequency (MHz);Magnitude [Arb Units]");
    // hm->SetStats(0);
    // hm->SetName("residualFFT");
    // hm->Scale(1.0 / hm->Integral());

    
    TCanvas *c1 = new TCanvas();
    c1->SetGrid();
    fftResidual->SetAxisRange(0,end/2,"X");
    fftResidual->SetTitle("Magnitude of the 1st transform");
    fftResidual->Draw();
    c1->SaveAs("figures/FFT_magtest.pdf");
    // hm->SetTitle("Magnitude of the 1st transform");
       // hm->Draw();
    //c1->SaveAs("figures/FFT_magtest.pdf");
    //NOTE: for "real" frequencies you have to divide the x-axes range with the range of your function
    //(in this case 4*Pi); y-axes has to be rescaled by a factor of 1/SQRT(n) to be right: this is not done automatically!
  
    //    hm->SetStats(kFALSE);
    // hm->GetXaxis()->SetLabelSize(0.05);
    // hm->GetYaxis()->SetLabelSize(0.05);
    // c1_3->cd();
    //Look at the phase of the output
    TH1 *hp = 0;
    hp = hsin->FFT(hp, "PH");

    TCanvas *c2 = new TCanvas();
    // gpad -> SetGrid();
    hp->SetTitle("Phase of the 1st transform");
    hp->Draw();
    // hp->SetStats(kFALSE);
    // hp->GetXaxis()->SetLabelSize(0.05);
    // hp->GetYaxis()->SetLabelSize(0.05);
    c2 ->SaveAs("figures/FFT_phtest.pdf");

 return 0;
}
