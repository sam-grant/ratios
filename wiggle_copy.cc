#include <iostream>
#include <fstream>
#include <iomanip>
//#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TFile.h"
#include "TLine.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TF1.h"
#include "TH1F.h"
#include "TF2.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TVirtualFFT.h"
#include "TFitter.h"
#include "TVirtualFitter.h"
#include "TFitResult.h"
#include <TMatrixDSym.h>
#include <cmath>
#include <math.h>
#include <iostream>
#include <vector>
#include <sstream>
#include <string>

using namespace std;

// Function for producing wiggle plot:
// All units are in ns

//Double_t totalTime = 30.0 * 4200.0; // Global

//Double_t binWidth = 150; // Time for muons to rotate once 

//Need to draw two points from pseudo experiments and substract

Double_t wiggle(Double_t *x, Double_t *par){
    
  Double_t time     = x[0];             // Leave time values free
  double tau       = 2.2;              // Lifetime of the muon at rest (ns)
  double gamma     = 29.3;              // Magic gamma
  double A         = 0.05;              // Amplitude
  double omega     = 2 * TMath::Pi() / 4.2;   // Time of a single wiggle ~ 4200 ns (omega=2pi/t)
  double omega2     = 2 * TMath::Pi() / 0.5; // Random additional high frequency with T = 10 ns
  double omega3     = 2 * TMath::Pi() / 10; // Random additional low frequency with T = 10000 ns
  double phase     = TMath::Pi()/2;     //Phase angle 
  double N = par[0];                        

  /* Decide whether to add frequencies together or not */
  //   Double_t Npositrons =  N * exp(- time / (tau * gamma) ) * (1 + A * cos ( (omega * time) + phase ));
  //  Double_t Npositrons =  N * exp(- time / (tau * gamma) ) * (1 + A * (cos ( (omega * time) + phase ) + cos ( (omega3 * time) + phase) ));//+ cos( (omega3 * time) + phase)));
  Double_t Npositrons =  N * (1 + A * (cos ( (omega * time) + phase ) + cos ( (omega3 * time) + phase) ));
  return Npositrons;
}

//Beckys analytic integral for the wiggle. 

Double_t integral(Double_t minTime, Double_t maxTime){

  double tau       = 2200;              // Lifetime of the muon at rest (ns)
  double gamma     = 29.3;              // Magic gamma
  double A         = 0.05;              // Amplitude
  double omega     = 2 * TMath::Pi() / 4200;   // Time of a single wiggle ~ 4200 ns (omega=2pi/t)
  double phase     = TMath::Pi()/2;
  double tau_gamma = tau * gamma;       //Making life easier
  double exp_min = exp(-minTime/tau_gamma);
  double exp_max = exp(-maxTime/tau_gamma);
  double arg_min = omega * minTime + phase;
  double arg_max = omega * maxTime + phase;
  double firstPart = tau_gamma * (exp_min - exp_max);
  double secondPart_const = A / ( (1/(tau_gamma*tau_gamma)) + omega*omega );
  double secondPart_min = exp_min * (omega * sin(arg_min) - (1/tau_gamma)*cos(arg_min));
  double secondPart_max = exp_max * (omega * sin(arg_max) - (1/tau_gamma)*cos(arg_max));
  double secondPart = secondPart_const * (secondPart_max - secondPart_min);

  double integral = firstPart + secondPart;

  return integral;
  
  }


// Need to insert a hi and low frequency 
int main() {

  const double binWidth = 0.15;
  const double totalTime = 30.0 * 4.2; // 30 wiggles

  TFile* file = new TFile("PseudoExp.root", "RECREATE"); // Create a ROOT file containing all pseudo experiments
  
  double Ntot = 2e+6;
  //Optional, set number of events. 
  // cout << "Please enter number of events: ";
  // cin >> Ntot;
  
  TF1 *f1 = new TF1("f1", wiggle, 0, totalTime, 1); // Book ideal wiggle plot
  f1 -> SetNpx(10000); // Prevent undersampling and make nice smooth plot SG
  f1->SetParameter(0,1); //  Force N to be 1. Else the function integral goes to zero. 
    
  double functionIntegral = f1->Integral(0,totalTime); 
  // double analyticIntegral = integral(0,totalTime);

  double N0 = Ntot/functionIntegral; // N0 the value of N at t = 0
  // double N0 = Ntot/analyticIntegral;
  cout << "N0 = " << N0 << endl;

  f1->SetParameter(0,N0); // Redfine normalisation to be N0, which normalises these 30 wiggles to the full 60 hr dataset
  
 // Loop over full time of experiment
  
  Double_t par[] = {1}; //What is this?
  Double_t Nfunc[int(totalTime / binWidth)]; // Array containing a number of elements equal to the number of rotations made. 
  double integralSum = 0;
  //  double analyticIntegralSum = 0;
  double functionIntegralSum = 0; 

  for (int ibin = 0; ibin < (totalTime / binWidth); ibin++){ //Loop over bins, or rotations. 

    double time = binWidth * ibin + binWidth/2;
    //Nfunc[ibin] = N0 * binWidth * wiggle(&time, par); // Time distribution for decay positrons

    //cout << "time = " << time << endl;
    
    /* The pseudo experiments draw from the integral not directly from the function.
     Confusing, when calculating binIntegral we use ROOTs integral, when calculating N0 we use the analytic integral. */


    double functionBinIntegral = f1->Integral(time-(binWidth/2), time+(binWidth/2));
    //    double analyticBinIntegral = integral(time-(binWidth/2), time+(binWidth/2));
    functionIntegralSum += functionBinIntegral;
    // analyticIntegralSum += analyticBinIntegral;
    //Nfunc[ibin] = analyticBinIntegral;// * N0;
    Nfunc[ibin] = functionBinIntegral; //why are we using the function integral now, works if there are bins??

    double centreValue = wiggle(&time, par);
    // f1->Write();
    
  }

    cout << "Function Integral = " << functionIntegral << ", Sum = " << functionIntegralSum << endl;
    //cout << "Integral over all bins = " << integralSum << endl;
    //cout << " *************** Analytic function integral = " << analyticIntegral << ", sum = " << analyticIntegralSum << endl;
    
  TRandom3 *r = new TRandom3(54321); // Set seed to get same random numbers (0 means no seed)

  // Loop over number of psuedo experiments (choosen by the user)

  int N_Pseudo_Exp;
  
  cout << "Please enter the number of pseudo experiements: " << endl;
  cin >>  N_Pseudo_Exp;   

  for (int i=0; i <= (N_Pseudo_Exp); i++){

    // Using stringstream to accommodate the number of pseudo experiements when naming and giving titles to plots
    
    string Pseudo_num;
    ostringstream conv;
    conv << i;
    string t0name  = "Pseudo0_" + conv.str();
    string t0title = "Pseudo0 Wiggle Plot " + conv.str();
    string t1name  = "PseudoA_" + conv.str();
    string t1title = "PseudoA Wiggle Plot " + conv.str();
    string t2name  = "PseudoB_" + conv.str();
    string t2title = "PseudoB Wiggle Plot " + conv.str();
    string t3name  = "PseudoC_" + conv.str();
    string t3title = "PseudoC Wiggle Plot " + conv.str();
    string t4name  = "PseudoD_" + conv.str();
    string t4title = "PseudoD Wiggle Plot " + conv.str();
    string tTotname  = "PseudoTot_" + conv.str();
    string tTottitle = "PseudoTot Wiggle Plot " + conv.str();
    
    string tDiffname  = "Pseudo_diff_" + conv.str();
    string tDifftitle = "Wiggle Difference Plot " + conv.str();

    // Book a histogram for pseudo wiggle plot
    TH1D* t0 = new TH1D("t0", "", int(totalTime / binWidth), 0, totalTime); // Undivided histogram for testing. 
    TH1D* t1 = new TH1D("t1", "", int(totalTime / binWidth), 0, totalTime);
    TH1D* t2 = new TH1D("t2", "", int(totalTime / binWidth), 0, totalTime);
    TH1D* t3 = new TH1D("t3", "", int(totalTime / binWidth), 0, totalTime);
    TH1D* t4 = new TH1D("t4", "", int(totalTime / binWidth), 0, totalTime);
    TH1D* tTot = new TH1D("tT0", "", int(totalTime / binWidth), 0, totalTime);

    t0 -> SetNameTitle(t0name.c_str(), t0title.c_str());
    t0 -> GetXaxis() -> SetTitle("Time (#us)");
    t0 -> GetYaxis() -> SetTitle("Number of Positrons N");

    t1 -> SetNameTitle(t1name.c_str(), t1title.c_str());
    t1 -> GetXaxis() -> SetTitle("Time (#us)");
    t1 -> GetYaxis() -> SetTitle("Number of Positrons N");

    t2 -> SetNameTitle(t2name.c_str(), t2title.c_str());
    t2 -> GetXaxis() -> SetTitle("Time (#us)");
    t2 -> GetYaxis() -> SetTitle("Number of Positrons N");

    t3 -> SetNameTitle(t3name.c_str(), t3title.c_str());
    t3 -> GetXaxis() -> SetTitle("Time (#us)");
    t3 -> GetYaxis() -> SetTitle("Number of Positrons N");

    t4 -> SetNameTitle(t4name.c_str(), t4title.c_str());
    t4 -> GetXaxis() -> SetTitle("Time (#us)");
    t4 -> GetYaxis() -> SetTitle("Number of Positrons N");

    tTot -> SetNameTitle(tTotname.c_str(), tTottitle.c_str());
    tTot -> GetXaxis() -> SetTitle("Time (#us)");
    tTot -> GetYaxis() -> SetTitle("Number of Positrons N");

    // Book a histogram for difference between ideal wiggle and pseudo wiggle plots
    
    TH1D* tDiff = new TH1D("", "", 100, -10, 10);
    
    tDiff -> SetNameTitle(tDiffname.c_str(), tDifftitle.c_str());
    tDiff -> GetXaxis() -> SetTitle("N/#sqrt{N}");
    tDiff -> GetYaxis() -> SetTitle("Number of Entries");    
    
    for (int ibin = 0; ibin < (totalTime / binWidth); ibin++){

      double nHits = Nfunc[ibin]/4.0; // 4 histograms, should be random, plot Nfunc
    
      double error = sqrt( Nfunc[ibin]/4.0 );

      Double_t PseudoWiggle0 =  Nfunc[ibin]; // Undivided histogram, no smearing. 
      Double_t PseudoWiggle1 = r->Gaus(nHits, error);  // Gaussian smearing
      Double_t PseudoWiggle2 = r->Gaus(nHits, error);  // Gaussian smearing
      Double_t PseudoWiggle3 = r->Gaus(nHits, error);  // Gaussian smearing
      Double_t PseudoWiggle4 = r->Gaus(nHits, error);  // Gaussian smearing
      
      Double_t PseudoWiggle = PseudoWiggle1 + PseudoWiggle2 + PseudoWiggle3 + PseudoWiggle4;

      t0 -> SetBinContent(ibin+1, PseudoWiggle0);
      t1 -> SetBinContent(ibin+1, PseudoWiggle1);
      t2 -> SetBinContent(ibin+1, PseudoWiggle2);
      t3 -> SetBinContent(ibin+1, PseudoWiggle3);
      t4 -> SetBinContent(ibin+1, PseudoWiggle4);
      
      tTot -> SetBinContent(ibin+1, PseudoWiggle);

      
      Double_t diff =  (Nfunc[ibin] - PseudoWiggle) / (sqrt(Nfunc[ibin])); // Difference between ideal and pseudo wiggle
      tDiff -> Fill(diff);
      
    }

    // Write TH1F to ROOT file
    // f1->Write();
    t0 -> Write();
    t1 -> Write();
    t2 -> Write();
    t3 -> Write();
    t4 -> Write();
    tTot -> Write();
    tDiff -> Write();

    // delete f1;
    delete t0;
    delete t1;
    delete t2;
    delete t3;
    delete t4;
    delete tTot;
    delete tDiff;

  }

  
  // Draw ideal wiggle plot
  
   f1->Write();
   TCanvas *c1 = new TCanvas();
  gPad->SetGrid();
  f1 -> SetTitle("Functional Wiggle Plot;Time [#us];Normalised Units");
  f1 -> Draw();
  c1 -> SaveAs("figures/IdealFit.pdf");

  delete c1;
  delete f1;
  file ->Close();
  
}
