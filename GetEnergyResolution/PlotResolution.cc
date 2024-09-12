#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TH1F.h"
#include "TPad.h"
#include "TRandom.h"
#include "TMath.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "LuxeStyle.h"
#include "LuxeStyle.C"
#include "LuxeLabels.h"
#include "LuxeLabels.C"

using namespace std;


string const PROCESSOR_NAME = "EcalEResolution";
string const FILE_LIST = "con_smearing.list";
string ecalType[3] = {"smeared", "pixelised", "monolithic"};
// string ecalType[3] = {"digitised", "pixelised", "monolithic"};

vector<string> GetFileNames(string const path) {
    ifstream filelist(path);
    string filename;
    vector<string> filenames;
    while (getline(filelist, filename)) {
        filenames.push_back(filename);
    }
    return filenames;
}

float GetRunEnergy(TFile* rootfile) {
    TH1F* runInfo = (TH1F*) rootfile->Get(Form("%s/_runInfo", PROCESSOR_NAME.c_str()));
    float energy = runInfo->GetBinContent(1);
    return energy;
}

vector<TF1*> GetFitFuncsHits(TFile* rootfile) {
    vector<TF1*> funcs;
    TH1D* hists[2]; TF1* func;
    hists[0] = (TH1D*) rootfile->Get(Form("%s/digit_evHitsHist", PROCESSOR_NAME.c_str()));
    hists[1] = (TH1D*) rootfile->Get(Form("%s/pixel_evHitsHist", PROCESSOR_NAME.c_str()));
    for (int i=0; i<2; i++) {
        func = (TF1*) hists[i]->GetListOfFunctions()->FindObject("gaus");
        funcs.push_back(func);
    }
    return funcs;
}

vector<TF1*> GetFitFuncsEnergy(TFile* rootfile) {
    vector<TF1*> funcs;
    TH1D* hists[3]; TF1* func;
    hists[0] = (TH1D*) rootfile->Get(Form("%s/digit_evEnergyHist", PROCESSOR_NAME.c_str()));
    hists[1] = (TH1D*) rootfile->Get(Form("%s/pixel_evEnergyHist", PROCESSOR_NAME.c_str()));
    hists[2] = (TH1D*) rootfile->Get(Form("%s/mono_evEnergyHist", PROCESSOR_NAME.c_str()));
    for (int i=0; i<3; i++) {
        func = (TF1*) hists[i]->GetListOfFunctions()->FindObject("gaus");
        funcs.push_back(func);
    }
    return funcs;
}

void PlotResolution() {
    SetLuxeStyle();
    
    uint plotMarker[3] = {5,2,4};// "x" "+" "o" markers
    map<float, double> energyMus[3], energySigmas[3];
    map<float, double> hitsMus[2], hitsSigmas[2];
    
    vector<string> filenames = GetFileNames(FILE_LIST);
    
    for (string filename : filenames) {
        cout<<"Processing file "<< filename <<endl;
        TFile* rootfile = TFile::Open(filename.c_str());
        float energy = GetRunEnergy(rootfile);
        
        vector<TF1*> energyfuncs = GetFitFuncsEnergy(rootfile);
        for (int i=0; i<3; i++) {
            energyMus[i][energy] = energyfuncs[i]->GetParameter(1);
            energySigmas[i][energy] = energyfuncs[i]->GetParameter(2);
        }
        vector<TF1*> hitsfuncs = GetFitFuncsHits(rootfile);
        for (int i=0; i<2; i++) {
            hitsMus[i][energy] = hitsfuncs[i]->GetParameter(1);
            hitsSigmas[i][energy] = hitsfuncs[i]->GetParameter(2);
        }
        rootfile->Close();
    }
    
    TCanvas* canvas = new TCanvas("canvas", "", 800, 600);
    // TCanvas* canvLin = new TCanvas("canvLin", "", 800, 600);
    TCanvas* canvRes = new TCanvas("canvRes", "", 800, 600);
    TCanvas* canvResiSqrt = new TCanvas("canvResiSqrt", "", 800, 600);
    TF1* funcLin = new TF1("funcLin", "1/[0]*x^[1]");
    TF1* funcRes = new TF1("funcRes", "sqrt([0]^2/x + ([1]/x)^2 + [2]^2)");
    TF1* funcResiSqrt = new TF1("funcResiSqrt", "sqrt(([0]*x)^2 + ([1]*x^2)^2 + [2]^2)");
    funcLin->SetParNames("p", "n");
    funcLin->SetParameters(1., 1.); funcLin->SetParLimits(1, 0.1, 10);
    funcLin->SetRange(1.0, 16.0); funcLin->SetLineStyle(9); funcLin->SetLineColor(kGray); funcLin->SetLineWidth(2);
    funcRes->SetParNames("a", "b", "c");
    funcRes->SetParameters(0.2, 0., 0.);
    funcRes->SetRange(1.0, 16.0); funcRes->SetLineStyle(9); funcRes->SetLineColor(kGray); funcRes->SetLineWidth(2);
    funcResiSqrt->SetParNames("a", "b", "c");
    funcResiSqrt->SetParameters(0.2, 0., 0.);
    funcResiSqrt->SetRange(0.2, 1.0); funcResiSqrt->SetLineStyle(9); funcResiSqrt->SetLineColor(kGray); funcResiSqrt->SetLineWidth(2);
    
    TGraph *energyRes[3], *energyResiSqrt[3]; TGraphErrors *energyLin[3];
    for (int i=0; i<3; i++) {
        energyLin[i] = new TGraphErrors();
        energyLin[i]->SetTitle(Form("ECAL-E %s; E_{0} [GeV]; E_{dep} [GeV]", ecalType[i].c_str()));
        energyRes[i] = new TGraph();
        energyRes[i]->SetTitle(Form("ECAL-E %s energy resolution; E_{0} [GeV]; #sigma_{E}/E", ecalType[i].c_str()));
        energyResiSqrt[i] = new TGraph();
        energyResiSqrt[i]->SetTitle(Form("ECAL-E %s energy resolution; 1/#sqrt{E_{0}} [GeV^{-1/2}]; #sigma_{E}/E", ecalType[i].c_str()));
        
        for (const auto&[energy, mu] : energyMus[i]) {
            energyLin[i]->AddPoint(energy, mu);
            energyLin[i]->SetPointError(energyLin[i]->GetN()-1, 0., energySigmas[i][energy]);
            energyRes[i]->AddPoint(energy, energySigmas[i][energy]/mu);
            energyResiSqrt[i]->AddPoint(1./TMath::Sqrt(energy), energySigmas[i][energy]/mu);
        }
        
        energyLin[i]->GetXaxis()->SetLimits(0, 16.5);
        energyLin[i]->SetMinimum(0); energyLin[i]->SetMaximum(0.175);
        if (i==0) {
            energyLin[i]->GetYaxis()->SetTitle("E_{dep} [MIP]");
            energyLin[i]->SetMinimum(0); energyLin[i]->SetMaximum(1250);
        }
        energyRes[i]->GetXaxis()->SetLimits(0, 16.5);
        energyRes[i]->SetMinimum(0); energyRes[i]->SetMaximum(0.25);
        energyResiSqrt[i]->GetXaxis()->SetLimits(0., 1.0);
        energyResiSqrt[i]->SetMinimum(0); energyResiSqrt[i]->SetMaximum(0.25);
        energyLin[i]->SetMarkerStyle(plotMarker[i]);
        energyRes[i]->SetMarkerStyle(plotMarker[i]);
        energyResiSqrt[i]->SetMarkerStyle(plotMarker[i]);
        
        if (i==0) {
            canvRes->cd();
            energyRes[i]->Draw("AP");
            canvResiSqrt->cd();
            energyResiSqrt[i]->Draw("AP");
        } else {
            canvRes->cd();
            energyRes[i]->Draw("P SAME");
            canvResiSqrt->cd();
            energyResiSqrt[i]->Draw("P SAME");
        }
    }
    // auto legendRes = new TLegend(0.2,0.2,0.45,0.45);
    // legendRes->SetHeader("Energy resolution"); // option "C" allows to center the header
    // legendRes->AddEntry(energyRes[2], Form("%s (energy)", ecalType[2].c_str()), "p");
    // legendRes->AddEntry(energyRes[1], Form("%s (energy)", ecalType[1].c_str()), "p");
    // legendRes->AddEntry(energyRes[0], Form("%s (energy)", ecalType[0].c_str()), "p");
    // canvRes->cd();
    // legendRes->Draw();
    // canvRes->Print("energy_resolution_inverse_sqrt_smeared.pdf");
    // canvRes->Clear();

    // auto legendResiSqrt = new TLegend(0.2,0.6,0.45,0.85);
    // legendResiSqrt->SetHeader("Energy resolution"); // option "C" allows to center the header
    // legendResiSqrt->AddEntry(energyResiSqrt[2], Form("%s (energy)", ecalType[2].c_str()), "p");
    // legendResiSqrt->AddEntry(energyResiSqrt[1], Form("%s (energy)", ecalType[1].c_str()), "p");
    // legendResiSqrt->AddEntry(energyResiSqrt[0], Form("%s (energy)", ecalType[0].c_str()), "p");
    // canvResiSqrt->cd();
    // legendResiSqrt->Draw();
    // canvResiSqrt->Print("energy_resolution_smeared.pdf");
    // canvResiSqrt->Clear();
    
    TGraph *hitsRes[2], *hitsResiSqrt[2]; TGraphErrors *hitsLin[2];
    for (int i=0; i<2; i++) {
        hitsLin[i] = new TGraphErrors();
        hitsLin[i]->SetTitle(Form("ECAL-E %s; E_{0} [GeV]; Hit", ecalType[i].c_str()));
        hitsRes[i] = new TGraph();
        hitsRes[i]->SetTitle(Form("ECAL-E %s hit resolution; E_{0} [GeV]; #sigma_{E}/E", ecalType[i].c_str()));
        hitsResiSqrt[i] = new TGraph();
        hitsResiSqrt[i]->SetTitle(Form("ECAL-E %s hit resolution; 1/#sqrt{E_{0}} [GeV^{-1/2}]; #sigma_{E}/E", ecalType[i].c_str()));
        for (const auto&[energy, mu] : hitsMus[i]) {
            hitsLin[i]->AddPoint(energy, mu);
            hitsLin[i]->SetPointError(hitsLin[i]->GetN()-1, 0., hitsSigmas[i][energy]);
            hitsRes[i]->AddPoint(energy, hitsSigmas[i][energy]/mu);
            hitsResiSqrt[i]->AddPoint(1./TMath::Sqrt(energy), hitsSigmas[i][energy]/mu);
        }
        hitsLin[i]->GetXaxis()->SetLimits(0, 16.5);
        hitsLin[i]->SetMinimum(0); hitsLin[i]->SetMaximum(250);
        hitsRes[i]->GetXaxis()->SetLimits(0, 16.5);
        hitsRes[i]->SetMinimum(0); hitsRes[i]->SetMaximum(0.25);
        hitsResiSqrt[i]->GetXaxis()->SetLimits(0., 1.0);
        hitsResiSqrt[i]->SetMinimum(0); hitsResiSqrt[i]->SetMaximum(0.25);
        hitsLin[i]->SetMarkerStyle(plotMarker[i]);      hitsLin[i]->SetMarkerColor(kRed);
        hitsRes[i]->SetMarkerStyle(plotMarker[i]);      hitsRes[i]->SetMarkerColor(kRed);
        hitsResiSqrt[i]->SetMarkerStyle(plotMarker[i]); hitsResiSqrt[i]->SetMarkerColor(kRed);
        
        // if (i==0) {
        //     canvRes->cd();
        //     hitsRes[i]->Draw("AP");
        //     canvResiSqrt->cd();
        //     hitsResiSqrt[i]->Draw("AP");
        // } else {
            canvRes->cd();
            hitsRes[i]->Draw("P SAME");
            canvResiSqrt->cd();
            hitsResiSqrt[i]->Draw("P SAME");
        // }
    }
    auto legendRes = new TLegend(0.2,0.2,0.45,0.45);
    legendRes->AddEntry(energyRes[2], Form("%s (energy)", ecalType[2].c_str()), "p");
    legendRes->AddEntry(energyRes[1], Form("%s (energy)", ecalType[1].c_str()), "p");
    legendRes->AddEntry(energyRes[0], Form("%s (energy)", ecalType[0].c_str()), "p");
    legendRes->AddEntry(hitsRes[1], Form("%s (hit)", ecalType[1].c_str()), "p");
    legendRes->AddEntry(hitsRes[0], Form("%s (hit)", ecalType[0].c_str()), "p");
    canvRes->cd();
    legendRes->Draw();
    canvRes->Print(Form("resolution_%s.pdf", ecalType[0].c_str()));
    canvRes->Clear();

    auto legendResiSqrt = new TLegend(0.2,0.6,0.45,0.85);
    legendResiSqrt->AddEntry(energyResiSqrt[2], Form("%s (energy)", ecalType[2].c_str()), "p");
    legendResiSqrt->AddEntry(energyResiSqrt[1], Form("%s (energy)", ecalType[1].c_str()), "p");
    legendResiSqrt->AddEntry(energyResiSqrt[0], Form("%s (energy)", ecalType[0].c_str()), "p");
    legendResiSqrt->AddEntry(hitsResiSqrt[1], Form("%s (hit)", ecalType[1].c_str()), "p");
    legendResiSqrt->AddEntry(hitsResiSqrt[0], Form("%s (hit)", ecalType[0].c_str()), "p");
    canvResiSqrt->cd();
    legendResiSqrt->Draw();
    canvResiSqrt->Print(Form("resolution_inverse_sqrt_%s.pdf", ecalType[0].c_str()));
    canvResiSqrt->Clear();

    // legendRes->SetHeader("Hit resolution"); // option "C" allows to center the header
    // canvRes->cd();
    // legendRes->Draw();
    // canvRes->Print("hits_resolution_inverse_sqrt.pdf");
    // canvRes->Clear();

    // legendResiSqrt->SetHeader("Hit resolution"); // option "C" allows to center the header
    // canvResiSqrt->cd();
    // legendResiSqrt->Draw();
    // canvResiSqrt->Print("hits_resolution.pdf");
    // canvResiSqrt->Clear();


    for (int i=0; i<3; i++) {
        canvas->cd();
        energyLin[i]->Draw("AP");
        energyLin[i]->Fit(funcLin, "R");
        gStyle->SetOptFit(1);//probability, chisquare, error, value
        canvas->Print(Form("energy_lin_ecale_%s.pdf", ecalType[i].c_str()));
        canvas->Clear();
        energyRes[i]->Draw("AP");
        energyRes[i]->Fit(funcRes, "R");
        gStyle->SetOptFit(1);//probability, chisquare, error, value
        canvas->Print(Form("energy_res_ecale_%s.pdf", ecalType[i].c_str()));
        canvas->Clear();
        energyResiSqrt[i]->Draw("AP");
        energyResiSqrt[i]->Fit(funcResiSqrt, "R");
        gStyle->SetOptFit(1);//probability, chisquare, error, value
        canvas->Print(Form("energy_res_inverse_sqrt_ecale_%s.pdf", ecalType[i].c_str()));
        canvas->Clear();
    }
    
    for (int i=0; i<2; i++) {
        canvas->cd();
        hitsLin[i]->Draw("AP");
        hitsLin[i]->Fit(funcLin, "R");
        gStyle->SetOptFit(1);//probability, chisquare, error, value
        canvas->Print(Form("hit_lin_ecale_%s.pdf", ecalType[i].c_str()));
        canvas->Clear();
        hitsRes[i]->Draw("AP");
        hitsRes[i]->Fit(funcRes, "R");
        gStyle->SetOptFit(1);//probability, chisquare, error, value
        canvas->Print(Form("hit_res_ecale_%s.pdf", ecalType[i].c_str()));
        canvas->Clear();
        hitsResiSqrt[i]->Draw("AP");
        hitsResiSqrt[i]->Fit(funcResiSqrt, "R");
        gStyle->SetOptFit(1);//probability, chisquare, error, value
        canvas->Print(Form("hit_res_inverse_sqrt_ecale_%s.pdf", ecalType[i].c_str()));
        canvas->Clear();
    }
    return;
}