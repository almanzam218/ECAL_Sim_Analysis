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


double GetFitParamsHits(int beamEnergyMeV,int param){
    TFile f(Form("/lustre/ific.uv.es/prj/gl/abehep.flc/LUXE/ECALe_SimAnalysis/electron/EnergyResolution/ECALe_luxe_v0_QGSP_BERT_e-_%dMeV_0.root",beamEnergyMeV));
    TH1F *hitsHist = (TH1F*)f.Get("MyGetEnergyResolutionProcessor/_evHitsHist");
    hitsHist->Fit("gaus");
    TF1 *fit = hitsHist->GetFunction("gaus");
    hitsHist->Fit("gaus","","",fit->GetParameter(1)-fit->GetParameter(2)*2,fit->GetParameter(1)+fit->GetParameter(2)*2);
    TF1 *fit1 = hitsHist->GetFunction("gaus");

    double fitParams[3] = {fit1->GetParameter(0),fit1->GetParameter(1),fit1->GetParameter(2)};
    f.Close();

    return fitParams[param];
}
double GetFitParamsEnergy(int beamEnergyMeV,int param){
    TFile f(Form("/lustre/ific.uv.es/prj/gl/abehep.flc/LUXE/ECALe_SimAnalysis/electron/EnergyResolution/ECALe_luxe_v0_QGSP_BERT_e-_%dMeV_0.root",beamEnergyMeV));
    TH1F *energyHist = (TH1F*)f.Get("MyGetEnergyResolutionProcessor/_evEnergyHist");
    energyHist->Fit("gaus");
    TF1 *fit = energyHist->GetFunction("gaus");
    energyHist->Fit("gaus","","",fit->GetParameter(1)-fit->GetParameter(2)*2,fit->GetParameter(1)+fit->GetParameter(2)*2);
    TF1 *fit1 = energyHist->GetFunction("gaus");

    double fitParams[3] = {fit1->GetParameter(0),fit1->GetParameter(1),fit1->GetParameter(2)};
    f.Close();

    return fitParams[param];
}

void PlotEnergyResolution(){

    TCanvas *c1 = new TCanvas("c1","Energy Resolution",1920,0,1920,1000);
    TCanvas *c2 = new TCanvas("c2","Energy Linearity",1920,0,1920,1000);
    TCanvas *c3 = new TCanvas("c3","Hit Resolution",1920,0,1920,1000);
    TCanvas *c4 = new TCanvas("c4","Hit Linearity",1920,0,1920,1000);
    

    double EnergyResolution[7];
    double meanEnergy[7];
    double HitResolution[7];
    double meanHits[7];
    double beamEnergy[7];


    int h=0;
    for (int i = 0; i < 7; i++)
    {
        EnergyResolution[h] = GetFitParamsEnergy(1500+2000*i,2)/GetFitParamsEnergy(1500+2000*i,1);
        meanHits[h] = GetFitParamsHits(1500+2000*i,1);
        HitResolution[h] = GetFitParamsHits(1500+2000*i,2)/GetFitParamsHits(1500+2000*i,1);
        meanEnergy[h] = GetFitParamsEnergy(1500+2000*i,1);
        beamEnergy[h] = (1500+2000*i)/1000;// in GeV
        h++;
    }
    auto gEnergyRes = new TGraph(7,beamEnergy,EnergyResolution);
    gEnergyRes->SetTitle("Energy Resolution");
    gEnergyRes->GetXaxis()->SetTitle("Beam Energy (MeV)");
    gEnergyRes->GetYaxis()->SetTitle("#frac{#sigma}{E}");
    
    auto gEnergyLin = new TGraph(7,beamEnergy,meanEnergy);
    gEnergyLin->SetTitle("Energy Linearity");
    gEnergyLin->GetXaxis()->SetTitle("Beam Energy (MeV)");
    gEnergyLin->GetYaxis()->SetTitle("mean E");
    
    auto gHitRes = new TGraph(7,beamEnergy,HitResolution);
    gHitRes->SetTitle("Hit Resolution");
    gHitRes->GetXaxis()->SetTitle("Beam Energy (MeV)");
    gHitRes->GetYaxis()->SetTitle("#frac{#sigma}{E}");
    
    auto gHitLin = new TGraph(7,beamEnergy,meanHits);
    gHitLin->SetTitle("Hit Linearity");
    gHitLin->GetXaxis()->SetTitle("Beam Energy (MeV)");
    gHitLin->GetYaxis()->SetTitle("mean hits");

    c1->cd();
    gEnergyRes->Draw("AC*");
    TF1 *f1 = new TF1("f1","sqrt(([0]^2/x)+([1]/x)^2+[2]^2)");
    f1->SetParameters(0,0,0);
    f1->SetParNames ("a","b","c");
    gEnergyRes->Fit(f1,"","",1,13);
    gStyle->SetOptFit(1);
    
    c2->cd();
    gEnergyLin->Draw("AC*");
    TF1 *f2 = new TF1("f2","[0]*x^[1]");
    f2->SetParameters(0,0);
    gEnergyLin->Fit(f2,"","",1,13);
    gStyle->SetOptFit(1);

    c1->Print("/lustre/ific.uv.es/prj/gl/abehep.flc/LUXE/ECALe_SimAnalysis/electron/EnergyResolution/EnergyResolutionPlot_20240821.png");
    c2->Print("/lustre/ific.uv.es/prj/gl/abehep.flc/LUXE/ECALe_SimAnalysis/electron/EnergyResolution/EnergyLinearityPlot_20240821.png");

    c3->cd();
    gHitRes->Draw("AC*");
    TF1 *f3 = new TF1("f3","sqrt(([0]^2/x)+([1]/x)^2+[2]^2)");
    f3->SetParameters(0,0,0);
    f3->SetParNames ("a","b","c");
    gHitRes->Fit(f3,"","",1,13);
    gStyle->SetOptFit(1);
    
    c4->cd();
    gHitLin->Draw("AC*");
    TF1 *f4 = new TF1("f4","[0]*x^[1]");
    f4->SetParameters(10,-10);
    gHitLin->Fit(f4,"","",1,13);
    gStyle->SetOptFit(1);

    c3->Print("/lustre/ific.uv.es/prj/gl/abehep.flc/LUXE/ECALe_SimAnalysis/electron/EnergyResolution/HitResolutionPlot_20240821.png");
    c4->Print("/lustre/ific.uv.es/prj/gl/abehep.flc/LUXE/ECALe_SimAnalysis/electron/EnergyResolution/HitLinearityPlot_20240821.png");

    return;
}