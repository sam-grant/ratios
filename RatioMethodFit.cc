#include <iostream>
#include <fstream>
#include <iomanip>
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TFile.h"
#include "TLine.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TF1.h"
#include "TF2.h"
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
#include <TDirectory.h>
using namespace std;

Double_t myfunction(Double_t *x, Double_t *par){

  Double_t time   = x[0];
  Double_t A     = par[0];
  Double_t omega = par[1];
  Double_t phase = par[2];
  Double_t delta = 2.87E-4; // From literature SG //par[3]; // calc don't fit

  Double_t Ratio =  par[0] * cos((par[1] * x[0]) + par[2]) + delta;

  return Ratio;
  
}

int main(){

  double w = 1500;
  double h = 1000;
  TCanvas* c2 = new TCanvas("c2", "c2", w, h);
  
  // Reading Root file and extracting histograms (choosen by user)
  
  int N_Pseudo_Exp;
  
  cout << "Please enter the number of histograms you want to use for the three parameter fit: " << endl;
  cin >>  N_Pseudo_Exp;
  
  TFile *file = TFile::Open("RatioMethodData.root", "READ");
  
  if (file == 0) {
    
    // if we cannot open the file, print an error message and return immediatly
    printf("Error: cannot open PseudoExp.root");
    return 1;
    
  }
  
  TFile* file3 = new TFile("ThreeParameterFit.root", "RECREATE");  // Create a ROOT file containg three parameter fits of the Ratio Method plots

  TH1D* h_amp = new TH1D("amp","", 1000, 0.049, 0.051);
  TH1D* h_omega = new TH1D("omega","", 1000, 0.0014955, 0.0014965);
  TH1D* h_phase = new TH1D("phase","", 1000, 0, 6.2);
  
  TH1D* h_prob = new TH1D("p-value","", 110, 0, 1.1);
  TH1D* h_chiSqFit = new TH1D("chiSqFit","", 1000, 0, 1500.);
  TH1D* h_chiSqPerNDFFit = new TH1D("chiSqPerNDFFit","", 100, 0, 5.);
  TH1D* h_NDFFit = new TH1D("ndfFit","", 1000, 0 -0.5, 1000. - 0.5);
  
  TH1D* h_chiSq = new TH1D("chiSq","", 1000, 0, 1500.);
  TH1D* h_chiSqPerNDF = new TH1D("chiSqPerNDF","", 100, 0, 5.);
  TH1D* h_NDF = new TH1D("ndf","", 1000, 0 - 0.5, 1000. - 0.5);

  TF1 *fit = new TF1 ("fit", myfunction, 20000, 100000, 4);   // Three parameter fit

  fit -> SetParameter(0, 1.00*0.05);
  fit -> SetParameter(1, 1.00*(2*M_PI/4200.0));
  fit -> SetParameter(2, 1.00*3*M_PI/2.0);
  //fit -> SetParameter(3, 2.87E-4);
    
  fit -> SetParLimits(0, 0.01, 0.10);
  fit -> SetParLimits(1, 0.0010, 0.0020);
  fit -> SetParLimits(2, 0.0, 2*M_PI);
  //fit -> SetParLimits(3, 1E-3, -1E-3);

  fit -> SetNpx(10000);
  fit -> SetNDF(4);
    
  fit -> SetParName(0, "A");
  fit -> SetParName(1, "Omega");
  fit -> SetParName(2, "Phi");
  fit -> SetParName(3, "Delta");
  
  fit -> SetLineColor(kRed);
  fit -> SetLineStyle(1);
  fit -> SetLineWidth(1);

  for (int i(0); i < (N_Pseudo_Exp); i++){                         // Loop over number of pseudo exp

    string Pseudo_Exp_;
    ostringstream conv1;
    conv1 << i;
    string pseudo_num = "Pseudo_Exp_" + conv1.str();
    
      
    //Create a main directory
    TDirectory *MainDirec = file3->mkdir(pseudo_num.c_str());
    MainDirec->cd();

    // Using stringstream to accommodate the number of pseudo experiements when naming and giving titles to plots
    
    string RatDist_Pseudo_num;
    ostringstream conv2;
    conv2 << i;
    string name = "RatDist_Pseudo_" + conv2.str();    
    
    TH1D* h1 = (TH1D*)file -> Get(name.c_str());                // Gets the ratio method plot
      
    double BinWidth = h1 -> GetBinWidth(1);                     // Gets bin width of ratio method plot

    double BinNum = h1 -> GetNbinsX();                          // Gets bin number of ratio method plot

    double Total_Time = h1->GetBinLowEdge(BinNum) + BinWidth;  //2100 + BinWidth * BinNum;                      // Find total time
    
    //h1-> GetXaxis() -> SetRangeUser(0,4400);
    
    h1 -> Fit(fit, "EMR+");
    
    gStyle->SetOptFit(1111);

    h1 -> Write();

    //ResDirec ->cd();
    
    ostringstream count;
    count << i;
    string ResName = "Residual_Pseudo_" + count.str();    
    TH1D* h2 = new TH1D(ResName.c_str(),"", BinNum , h1->GetBinLowEdge(1),  h1->GetBinLowEdge(1) + Total_Time);

    ostringstream count2;
    count2 << i;
    string ResPullName = "Residual_Pull_Pseudo_" + count2.str();    
    TH1D* h3 = new TH1D(ResPullName.c_str(),"", BinNum , h1->GetBinLowEdge(1),  h1->GetBinLowEdge(1) + Total_Time);

    ostringstream count3;
    count3 << i;
    string ResPullDistName = "Residual_Pull_Distribution_Pseudo_" + count3.str();    
    TH1D* h4 = new TH1D(ResPullDistName.c_str(),"", 64, -15, 15);

    double chiSq = 0.0;
    int binsInChiSq = 0;
    
    //get our histogram and loop over bins
    for (int ibin(1); ibin < BinNum; ibin++){
      double time = h1->GetBinCenter(ibin);
      double measured = h1->GetBinContent(ibin);
      double error = h1->GetBinError(ibin);
      double fitted = fit->Eval(time);
      
      //cout << "time: " << time << ", meas: " << measured << ", err: " << error << ", fitted: " << fitted << "\n";
      
      // meas - fitted / error

      if (error > 0){
	h2 -> SetBinContent(ibin,(measured-fitted)); //residual SG

	h3 -> SetBinContent(ibin,(measured-fitted)/error); //pull SG

	//h2 -> SetBinContent(ibin,(measured-fitted)/error);
	//h3 -> SetBinContent(ibin,pow((measured-fitted),2)/pow(error,2));

	h4 -> Fill((measured-fitted)/error);
	//h4 -> Fill( pow((measured-fitted),2)/pow(error,2) );
	chiSq += (measured-fitted) / error;
	binsInChiSq++;
      }
      
    }

    
    h2 -> Write();
    h3 -> Write();
    h4 -> Write();

    h2->Delete();
    h3->Delete();
    h4->Delete();
    
    h_chiSq->Fill(chiSq);
    h_chiSqPerNDF->Fill(chiSq / double(binsInChiSq));
    h_NDF -> Fill(binsInChiSq);

    //from fit
    h_prob->Fill(fit->GetProb());
    h_chiSqFit->Fill(fit->GetChisquare());
    h_chiSqPerNDFFit->Fill( fit->GetChisquare() / double(fit->GetNDF()));
    h_NDFFit->Fill(fit->GetNDF());

    h_amp->Fill(fit->GetParameter(0));
    h_omega->Fill(fit->GetParameter(1));
    h_phase->Fill(fit->GetParameter(2));
    
    cout << "chi sq test, us: " << chiSq << " perndf: " << chiSq/ double(binsInChiSq) << " from fit: " << fit->GetChisquare() << ", per ndf " << fit->GetChisquare() / fit->GetNDF() << " with prob " << fit->GetProb() << ", total time: " << Total_Time << "\n"; 
    
    MainDirec->Close();
    
    // Get fitter
    
    //TVirtualFitter* f = TVirtualFitter::GetFitter();

    // Get number of total parameters in the fit
    
    //const int n = f -> GetNumberTotalParameters();
    
    //if (0 == f->GetCovarianceMatrix()) {
    //std::cout << "Fit failed, we therefore have no covariance matrix." <<
    //std::endl;
    //return 1;
    //}
    
    // Get covariance matrix
    
    //TMatrixDSym cov(n, f->GetCovarianceMatrix());
    
    // Print correlation matrix (conversion from covariance matrix to correlation matrix is done on the fly)
    
    //std::cout << "Correlation Matrix of the fit is: " << std::endl;
    
    //std::cout << setw(15) << right << " A " << setw(15) << right << " Omega " << setw(15) << right << " Phi " << std::endl;
    
    //for (int i = 0; i < n; ++i) {
      
    //switch(i){
	
    //case 0:
    //cout << setw(10) << left << "A" ;
    //break;
    //case 1:
    //cout << setw(10) << left << "Omega" ;
    //break;
    //case 2:
    //cout << setw(10) << left << "Phi" ;
    //break;
	
    //}
      
    //for (int j = 0; j < n; ++j)
    //std::cout << " " <<
    //std::setw(12) << std::setprecision(6) << 
    //cov(i, j) / std::sqrt(cov(i, i) * cov(j, j));
    //std::cout << std::endl;
      
    //}
    
    // Ideal values
    
    //double ideal_A = 0.05;
    //double ideal_omega = 2 * M_PI / 4200;
    //double ideal_phase = (M_PI/2) + (M_PI);
    
    // Fitted parameters

    //double fitted_A      = fit -> GetParameter(0);
    //double fitted_omega  = fit -> GetParameter(1);
    //double fitted_phase  = fit -> GetParameter(2);

    // Correlation plots and fits for fitted parameters

    //CorrFit1 = new TF1 ("CorrFit1", "pol1", 0.04999, 0.05001);   // A & Omega correlation fit
    //CorrFit1 -> SetParName(0, "C");
    //CorrFit1 -> SetParName(1, "Grad");

    //CorrFit2 = new TF1 ("CorrFit2", "pol1", 0.04999, 0.05001);   // A & Phi correlation fit
    //CorrFit2 -> SetParName(0, "C");
    //CorrFit2 -> SetParName(1, "Grad");

    //CorrFit3 = new TF1 ("CorrFit3", "pol1",  0.00149598, 0.00149601);   // Omega & Phi correlation fit
    //CorrFit3 -> SetParName(0, "C");
    //CorrFit3 -> SetParName(1, "Grad");
    
    //k1 -> Fill(fitted_A, fitted_omega);
    //k1 -> Fit("CorrFit1","R");
    //k1 -> SetTitle("Correlation between A_{Fitted} and #omega_{Fitted} (3 Parameter Fit)");
    //k1 -> GetXaxis() -> SetTitle("A_{Fitted}");
    //k1 -> GetYaxis() -> SetTitle("#omega_{Fitted}");
    //k1 -> SetMarkerStyle(20);
    //k1 -> SetMarkerSize(0.5);
    //k1 -> SetMarkerColor(kBlack);
    
    //k2 -> Fill(fitted_A, fitted_phase);
    //k2 -> Fit("CorrFit2","R");
    //k2 -> SetTitle("Correlation between A_{Fitted} and #phi_{Fitted} (3 Parameter Fit)");
    //k2 -> GetXaxis() -> SetTitle("A_{Fitted}");
    //k2 -> GetYaxis() -> SetTitle("#phi_{Fitted}");
    //k2 -> SetMarkerStyle(20);
    //k2 -> SetMarkerSize(0.5);
    //k2 -> SetMarkerColor(kBlack);
    
    //k3 -> Fill(fitted_omega, fitted_phase);
    //k3 -> Fit("CorrFit3","R");
    //k3 -> SetTitle("Correlation between #omega_{Fitted} and #phi_{Fitted} (3 Parameter Fit)");
    //k3 -> GetXaxis() -> SetTitle("#omega_{Fitted}");
    //k3 -> GetYaxis() -> SetTitle("#phi_{Fitted}");
    //k3 -> SetMarkerStyle(20);
    //k3 -> SetMarkerSize(0.5);
    //k3 -> SetMarkerColor(kBlack);

    // ChiSquare distribution histogram

    //double ChiSq = fit -> GetChisquare()/fit->GetNDF();
    //cs -> Fill(ChiSq);
    //cs -> SetTitle("#chi^{2} distribution (3 Parameter Fit)");
    //cs -> GetXaxis() -> SetTitle("#chi^{2}");
    //cs -> GetYaxis() -> SetTitle("Number of Entries");
    
    // Total difference histograms

    //double diff_A = (ideal_A - fitted_A);
    
    //g1 -> Fill(diff_A);
    //g1 -> SetTitle("Total difference in A (3 Parameter Fit)");
    //g1 -> GetXaxis() -> SetTitle(" A_{Ideal} - A_{Fitted} ");
    //g1 -> GetYaxis() -> SetTitle("Number of Entries");
    
    //double diff_omega = (ideal_omega - fitted_omega);
    
    //g2 -> Fill(diff_omega);
    //g2 -> SetTitle("Total difference in #omega (3 Parameter Fit)");
    //g2 -> GetXaxis() -> SetTitle("#omega_{Ideal} - #omega_{Fitted}");
    //g2 -> GetYaxis() -> SetTitle("Number of Entries");
      
    //double diff_phase = (ideal_phase - fitted_phase);
    
    //g3 -> Fill(diff_phase);
    //g3 -> SetTitle("Total difference in #phi (3 Parameter Fit)");
    //g3 -> GetXaxis() -> SetTitle("#phi_{Ideal} - #phi_{Fitted}");
    //g3 -> GetYaxis() -> SetTitle("Number of Entries");

    // Error on fitted parameters histograms

    //double ErrorPar0 = fit -> GetParError(0);
    
    //e1 -> Fill(ErrorPar0);
    //e1 -> SetTitle("Error on A_{Fitted} (3 Parameter Fit)");
    //e1 -> GetXaxis() -> SetTitle("Error");
    //e1 -> GetYaxis() -> SetTitle("Number of Entries");

    //double ErrorPar1 = fit -> GetParError(1);
    
    //e2 -> Fill(ErrorPar1);
    //e2 -> SetTitle("Error on #omega_{Fitted} (3 Parameter Fit)");
    //e2 -> GetXaxis() -> SetTitle("Error");
    //e2 -> GetYaxis() -> SetTitle("Number of Entries");

    //double ErrorPar2 = fit -> GetParError(2);
    
    //e3 -> Fill(ErrorPar2);
    //e3 -> SetTitle("Error on #phi_{Fitted} (3 Parameter Fit)");
    //e3 -> GetXaxis() -> SetTitle("Error");
    //e3 -> GetYaxis() -> SetTitle("Number of Entries");

    //    h1 -> Write();
    //h2 -> Write();
  }

  file3->cd();
  h_chiSq->Write();
  h_chiSqPerNDF->Write();
  h_NDF -> Write();
  h_prob->Write();
  h_chiSqFit->Write();
  h_chiSqPerNDFFit->Write();
  h_NDFFit->Write();
  h_amp->Write();
  h_omega->Write();
  h_phase->Write();

  h_chiSq->Delete();
  h_chiSqPerNDF->Delete();
  h_NDF->Delete();
  h_prob->Delete();
  h_chiSqFit->Delete();
  h_chiSqPerNDFFit->Delete();
  h_NDFFit->Delete();
  h_amp->Delete();
  h_omega->Delete();
  h_phase->Delete();

  //k1 -> Draw();
  //CorrFit1 -> Draw("SAME");
  //c2 -> SaveAs("correlation_A_omega_3PF.eps");
  
  //k2 -> Draw();
  //CorrFit2 -> Draw("SAME");
  //c2 -> SaveAs("correlation_A_phase_3PF.eps");
  
  //k3 -> Draw();
  //CorrFit3 -> Draw("SAME");
  //c2 -> SaveAs("correlation_omega_phase_3PF.eps");

  ///

  //cs -> Draw();
  //c2 -> SaveAs("ChiSquare_3PF.eps");
  
  ///
  
  //g1 -> Draw();
  //c2 -> SaveAs("Diff_A_3PF.eps");

  //g2 -> Draw();
  //c2 -> SaveAs("Diff_Omega_3PF.eps");
  
  //g3 -> Draw();
  //c2 -> SaveAs("Diff_Phase_3PF.eps");

  ///
  
  //e1 -> Draw();
  //c2 -> SaveAs("Error_A_3PF.eps");
  
  //e2 -> Draw();
  //c2 -> SaveAs("Error_Omega_3PF.eps");

  //e3 -> Draw();
  //c2 -> SaveAs("Error_Phase_3PF.eps");
  
}
