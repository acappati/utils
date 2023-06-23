// **********
// run with root -b -q minBiasXsec.C++
// **********

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "TGraphErrors.h"
#include "TF1.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TPaveText.h"


float lin_func(float x, float par0, float par1){

  return par0 + par1 * x;

}


void minBiasXsec(){

  const int n = 3;
  float energy[n]      = {7., 8., 13.}; // [TeV]
  float xsecMinBias[n] = {68.000, 69.400, 69.200}; // [mb]
  float en_run3 = 13.6; // [TeV]
  float x_min = 5.;
  float x_max = 15.;
  

  TCanvas* c = new TCanvas("c","",800,600);
  c->cd();
  TGraph* g = new TGraph(n, energy, xsecMinBias);
  g->SetMarkerStyle(20);
  g->SetTitle("");
  g->GetXaxis()->SetTitle("E [TeV]");
  g->GetYaxis()->SetTitle("xsec minBias [mb]");
  g->GetYaxis()->SetRangeUser(60.,75.);

  TF1 *f1 = new TF1("f1", "pol1", x_min, x_max);
  g->Fit("f1", "R"); // linear fit
  gStyle->SetOptFit(1111);
  g->Draw("APE");

  // eval function for run3 energy
  TF1 *fit = g->GetFunction("f1");
  float p0 = fit->GetParameter(0);
  float p1 = fit->GetParameter(1);

  float xsecMinBias_run3 = lin_func(en_run3, p0, p1);
  cout<<xsecMinBias_run3<<endl;

  c->SaveAs("minBiasXsec.png");
  
}
