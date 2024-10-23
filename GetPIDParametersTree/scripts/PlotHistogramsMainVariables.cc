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
void PlotHistogramsMainVariables(){

 
////GAMMA
   TFile fgamma("/lustre/ific.uv.es/prj/gl/abehep.flc/LUXE/ECALe_SimAnalysis/PIDParametersTrees/20240827_v1/histogramsIncluded/PIDParams_PixelDigiCluster_ECALe_luxe_v1_QGSP_BERT_gamma_0.5to10GeV.root");
   TTree *tgamma = (TTree*)fgamma.Get("ntp");
   Float_t MIP_Likeness, bar_z, hits_max_distance, sume_layer_5, mol,nhits;

    tgamma->SetBranchAddress("MIP_Likeness",&MIP_Likeness);
    tgamma->SetBranchAddress("bar_z",&bar_z);
    tgamma->SetBranchAddress("hits_max_distance",&hits_max_distance);
    tgamma->SetBranchAddress("sume_layer_5",&sume_layer_5);
    tgamma->SetBranchAddress("mol",&mol);
    tgamma->SetBranchAddress("nhit",&nhits);


   TH1F *MIP_LikenessHist;
   TH1F *bar_zHist;
   TH1D *nhitsHist;
   TH1F *hits_max_distanceHist;
   TH1F *sume_layer_5Hist;
   TH1F *molHist;

    MIP_LikenessHist = new TH1F("MIP_Likeness","MIP Likeness",50,0,1);
    bar_zHist = new TH1F("bar_z","barycenter z", 100, 0, 200);
    nhitsHist = new TH1D("nhits","Total hits", 210, -0.5, 209.5);
    hits_max_distanceHist = new TH1F("hits_max_distance","Hits Maximum Distance", 200 ,0 , 350);
    sume_layer_5Hist = new TH1F("sume_layer","Sum energy layer 5", 200, 0, 900);
    molHist = new TH1F("mol","Moliere radius",100,0, 150);
        
    for(int i=0;i<tgamma->GetEntries();i++){
        tgamma->GetEntry(i);
        MIP_LikenessHist->Fill(MIP_Likeness);
        bar_zHist->Fill(bar_z);
        nhitsHist->Fill(nhits);
        hits_max_distanceHist->Fill(hits_max_distance);
        sume_layer_5Hist->Fill(sume_layer_5);
        molHist->Fill(mol);   
    }
///NEUTRON
   TFile fneutron("/lustre/ific.uv.es/prj/gl/abehep.flc/LUXE/ECALe_SimAnalysis/PIDParametersTrees/20240827_v1/histogramsIncluded/PIDParams_PixelDigiCluster_ECALe_luxe_v1_QGSP_BERT_neutron_0.5to10GeV.root");
   TTree *tneutron = (TTree*)fneutron.Get("ntp");
   Float_t MIP_Likenessn, bar_zn, hits_max_distancen, sume_layer_5n, moln,nhitsn;

    tneutron->SetBranchAddress("MIP_Likeness",&MIP_Likenessn);
    tneutron->SetBranchAddress("bar_z",&bar_zn);
    tneutron->SetBranchAddress("hits_max_distance",&hits_max_distancen);
    tneutron->SetBranchAddress("sume_layer_5",&sume_layer_5n);
    tneutron->SetBranchAddress("mol",&moln);
    tneutron->SetBranchAddress("nhit",&nhitsn);


   TH1F *MIP_LikenessHistn;
   TH1F *bar_zHistn;
   TH1D *nhitsHistn;
   TH1F *hits_max_distanceHistn;
   TH1F *sume_layer_5Histn;
   TH1F *molHistn;

    MIP_LikenessHistn = new TH1F("MIP_Likenessn","MIP Likeness",50,0,1);
    bar_zHistn = new TH1F("bar_zn","barycenter z", 100, 0, 200);
    nhitsHistn = new TH1D("nhitsn","Total hits", 210, -0.5, 209.5);
    hits_max_distanceHistn = new TH1F("hits_max_distancen","Hits Maximum Distance", 200 ,0 , 350);
    sume_layer_5Histn = new TH1F("sume_layern","Sum energy layer 5", 200, 0, 900);
    molHistn = new TH1F("moln","Moliere radius",100,0, 150);
        
    for(int i=0;i<tneutron->GetEntries();i++){
        tneutron->GetEntry(i);
        MIP_LikenessHistn->Fill(MIP_Likenessn);
        bar_zHistn->Fill(bar_zn);
        nhitsHistn->Fill(nhitsn);
        hits_max_distanceHistn->Fill(hits_max_distancen);
        sume_layer_5Histn->Fill(sume_layer_5n);
        molHistn->Fill(moln);   
    }



///PI-
   TFile fpi("/lustre/ific.uv.es/prj/gl/abehep.flc/LUXE/ECALe_SimAnalysis/PIDParametersTrees/20240827_v1/histogramsIncluded/PIDParams_PixelDigiCluster_ECALe_luxe_v1_QGSP_BERT_pi-_0.5to10GeV.root");
   TTree *tpi = (TTree*)fpi.Get("ntp");
   Float_t MIP_Likenessp, bar_zp, hits_max_distancep, sume_layer_5p, molp,nhitsp;

    tpi->SetBranchAddress("MIP_Likeness",&MIP_Likenessp);
    tpi->SetBranchAddress("bar_z",&bar_zp);
    tpi->SetBranchAddress("hits_max_distance",&hits_max_distancep);
    tpi->SetBranchAddress("sume_layer_5",&sume_layer_5p);
    tpi->SetBranchAddress("mol",&molp);
    tpi->SetBranchAddress("nhit",&nhitsp);


   TH1F *MIP_LikenessHistp;
   TH1F *bar_zHistp;
   TH1D *nhitsHistp;
   TH1F *hits_max_distanceHistp;
   TH1F *sume_layer_5Histp;
   TH1F *molHistp;

    MIP_LikenessHistp = new TH1F("MIP_Likenessp","MIP Likeness",50,0,1);
    bar_zHistp = new TH1F("bar_zp","barycenter z", 100, 0, 200);
    nhitsHistp = new TH1D("nhitsp","Total hits", 210, -0.5, 209.5);
    hits_max_distanceHistp = new TH1F("hits_max_distancep","Hits Maximum Distance", 200 ,0 , 350);
    sume_layer_5Histp = new TH1F("sume_layerp","Sum energy layer 5", 200, 0, 900);
    molHistp = new TH1F("molp","Moliere radius",100,0, 150);
        
    for(int i=0;i<tpi->GetEntries();i++){
        tpi->GetEntry(i);
        MIP_LikenessHistp->Fill(MIP_Likenessp);
        bar_zHistp->Fill(bar_zp);
        nhitsHistp->Fill(nhitsp);
        hits_max_distanceHistp->Fill(hits_max_distancep);
        sume_layer_5Histp->Fill(sume_layer_5p);
        molHistp->Fill(molp);   
    }

    TCanvas *c1 = new TCanvas("c1","MIP Likeness",1920,0,1920,1000);
    TCanvas *c2 = new TCanvas("c2","Barycenter z",1920,0,1920,1000);
    TCanvas *c3 = new TCanvas("c3","Total Hits",1920,0,1920,1000);
    TCanvas *c4 = new TCanvas("c4","Hits maximum distance",1920,0,1920,1000);
    TCanvas *c5 = new TCanvas("c5","Sum energy layer 5",1920,0,1920,1000);
    TCanvas *c6 = new TCanvas("c6","Moliere radius",1920,0,1920,1000);

    c1->cd();
    MIP_LikenessHist->Draw();
    MIP_LikenessHistn->Draw("SAME");
    MIP_LikenessHistp->Draw("SAME");
    c1->Print("/lustre/ific.uv.es/prj/gl/abehep.flc/LUXE/ECALe_SimAnalysis/PIDParametersTrees/20240827_v1/histogramsIncluded/MIP_Likeness.png");
    c1->Print("/lustre/ific.uv.es/prj/gl/abehep.flc/LUXE/ECALe_SimAnalysis/PIDParametersTrees/20240827_v1/histogramsIncluded/MIP_Likeness.C");

    //c1->Clear();
    //c1->Close();
    fgamma.Close();
    fneutron.Close();
    fpi.Close();
    return;
}
