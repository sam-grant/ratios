/* Somewhere the test TVirtualFFT - SG 2018 */

#include <iostream>
#include <TVirtualFFT.h>
#include <TF1.h>
#include <TH1D.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TDirectory.h>
#include <string>

using namespace std;

int main() {

  cout << "Please enter the number of histograms you want to use for the FFT: " << endl;
  cin >>  N_Pseudo_Exp;

  TFile *file1 = TFile::Open("ThreeParameterFit.root", "READ");

  if (file == 0) {

    // if we cannot open the file, print an error message and return immediatly                                  
    printf("Error: cannot open root file");
    return 1;

  }

  TFile *file2 = TFile::Open("FFT.root", "RECREATE");

  
  for (int i(0); i < (N_Pseudo_Exp); i++){
  
     string Pseudo_Exp_;
    ostringstream conv1;
    conv1 << i;
    string pseudo_exp_num = "Pseudo_Exp_" + conv1.str();
    
    string Residual_Pseudo_;
    ostringstream conv2;
    conv2 << i;
    string name = "Residual_Pseudo_" + conv2.str();

    string FFT_Residual_;
    ostringstream conv3;
    conv3 << i;
    sring fft_name = "FFT_Residual_" + conv3.str();

    TH1D *h1 = (TH1D*)file1 -> Get(name.c_str());

  }

  
 return 0;
}
