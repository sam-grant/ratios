#include <iostream>
#include <TF1.h>
#include <TMath.h>

using namespace std;


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

  Double_t Npositrons = N * exp( - time / (tau * gamma)) * (1 + A * cos ( (omega * time) + phase)) ;
  //  Double_t Npositrons =  N * exp(- time / (tau * gamma) ) * (1 + A * (cos ( (omega * time) + phase ) + cos ( (omega2 * time) + phase) + cos( (omega3 * time) + phase)));
 return Npositrons;
}


int main() {

  const double totalTime = 30.0 * 4200.0;

  TF1 *f1 = new TF1("f1", Wiggle , 0 , totalTime , 1);  
  f1 -> SetNpx(10000);
  f1 -> SetParameter(0,1); // You have to force the normalisation to be one.

  double functionIntegral = f1 -> Integral(0, totalTime);
  cout<<functionIntegral<<endl;

  delete f1;
   return 0;

}
