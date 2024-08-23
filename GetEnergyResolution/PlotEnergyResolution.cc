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
    TFile f(Form("/lustre/ific.uv.es/prj/gl/abehep.flc/LUXE/ECALe_SimAnalysis/electron/EnergyResolution/ECALe_luxe_v0_QGSP_BERT_e-_%dMeV.root",beamEnergyMeV));
    TH1F *energyHist = (TH1F*)f.Get("MyGetEnergyResolutionProcessor/_evEnergyHist");
    energyHist->Fit("gaus");
    TF1 *fit = energyHist->GetFunction("gaus");

    double EnergyResolution = fit->GetParameter(2)/fit->GetParameter(1);
    f.Close();

    return EnergyResolution;
}

void PlotEnergyResolution(){

    TCanvas *c1 = new TCanvas("c1","Energy Resolution",1920,0,1920,1000);
    

    double Resolution[4];
    double beamEnergy[4];


    int h=0;
    for (int i = 0; i < 4; i++)
    {
        Resolution[h] = GetFitParams(1500+500*i);
        beamEnergy[h] = (1500+500*i);// in GeV
        h++;
    }
    auto gEnergyRes = new TGraph(4,beamEnergy,Resolution);

    c1->cd();
    gEnergyRes->Draw("AC*");
    c1->Print("testPlot.png");


    return;
}