#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TH1F.h"
#include "TPad.h"
#include "TRandom.h"
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>


double GetFitParams(int beamEnergyMeV){
    TFile f(Form("/lustre/ific.uv.es/prj/gl/abehep.flc/LUXE/ECALe_SimAnalysis/electron/EnergyResolution/ECALe_luxe_v0_QGSP_BERT_e-_%dMeV_0.root",beamEnergyMeV));
    TH1F *energyHist = (TH1F*)f.Get("MyGetEnergyResolutionProcessor/_evEnergyHist");
    energyHist->Fit("gaus");
    TF1 *fit = energyHist->GetFunction("gaus");

    double EnergyResolution = fit->GetParameter(2)/fit->GetParameter(1);
    f.Close();

    return EnergyResolution;
}

void PlotEnergyResolution(){

    TCanvas *c1 = new TCanvas("c1","Energy Resolution",1920,0,1920,1000);
    

    double Resolution[7];
    double beamEnergy[7];


    int h=0;
    for (int i = 0; i < 7; i++)
    {
        Resolution[h] = GetFitParams(1500+2000*i);
        beamEnergy[h] = (1500+2000*i);// in GeV
        h++;
    }
    auto gEnergyRes = new TGraph(7,beamEnergy,Resolution);
    gEnergyRes->SetTitle("Energy Resolution");
    gEnergyRes->GetXaxis()->SetTitle("Beam Energy (MeV)");
    gEnergyRes->GetYaxis()->SetTitle("#frac{#sigma}{E}");
    

    c1->cd();
    gEnergyRes->Draw("AC*");
    TF1 *f1 = new TF1("f1","sqrt(([0]^2/x)+([1]/x)^2+[2]^2)",1500,13500);
    f1->SetParameters(1,1,0);
    f1->SetParNames ("a","b","c");
    gEnergyRes->Fit(f1,"","",1500,13500);
    gStyle->SetOptFit(1);

    c1->Print("/lustre/ific.uv.es/prj/gl/abehep.flc/LUXE/ECALe_SimAnalysis/electron/EnergyResolution/EnergyResolutionPlot_20240821.png");


    return;
}