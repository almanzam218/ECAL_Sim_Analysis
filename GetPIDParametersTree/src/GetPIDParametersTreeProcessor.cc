#include "GetPIDParametersTreeProcessor.hh"
// #include "langaus.C"

// ROOT
#include "TStyle.h"
#include "TDirectory.h"
#include "TObjArray.h"
#include "TMath.h"

#include <math.h>
#include <iostream>
#include <fstream>
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TVectorD.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include <iostream>
#include <numeric>
// ----- include for verbosity dependent logging ---------
// #include "marlin/VerbosityLevels.h"
// #include "marlin/StringParameters.h"
// #define SLM streamlog_out(MESSAGE)

// using namespace std;
using namespace lcio ;
using namespace marlin ;
using EVENT::LCCollection;
using EVENT::MCParticle;
using EVENT::ReconstructedParticle;
using EVENT::Track;
using EVENT::Vertex;
using IMPL::LCRelationImpl;
using IMPL::ReconstructedParticleImpl;
using IMPL::TrackImpl;
using IMPL::TrackStateImpl;
using std::string;
using std::vector;
using UTIL::LCRelationNavigator;

GetPIDParametersTreeProcessor aGetPIDParametersTreeProcessor;

GetPIDParametersTreeProcessor::GetPIDParametersTreeProcessor() : Processor("GetPIDParametersTreeProcessor") {

	// modify processor description
	_description = "";

	// input collections
	registerInputCollection(LCIO::MCPARTICLE,"MCCollectionName",
                            "Name of the MC collection",
							_MCColName,
							// std::string("MCParticle"));
							std::string("Not configured in xml file"));
	registerInputCollection(LCIO::SIMCALORIMETERHIT, "ECALCollectionName",
							"Monolithic ECAL Hits Collection",
							_ECALColName,
							// std::string("SiEcalCollection"));
							std::string("Not configured in xml file"));
	registerInputCollection(LCIO::SIMCALORIMETERHIT, "PixelisedECALCollectionName",
							"Pixelised ECAL Hits Collection",
							_pECALColName,
							// std::string("PixelSiEcalCollection"));
							std::string("Not configured in xml file"));
	registerInputCollection(LCIO::CALORIMETERHIT, "DigitisedECALCollectionName",
							"Digitised ECAL Hits Collection",
							_dECALColName,
							// std::string("EcalCollection"));
							std::string("Not configured in xml file"));
}

GetPIDParametersTreeProcessor::~GetPIDParametersTreeProcessor() {}

void GetPIDParametersTreeProcessor::get_res(int &nhit, float &sume, float &weight, vector<float> * hit_energy, vector<int> *hit_slab, TVectorD W_thicknesses, vector<int> *hit_isMasked, bool &masked) {
    if(hit_energy->size() > 0){    
        //cout<<W_thicknesses.Min()<<endl;
        // First option: use the minimum of the used
        // Second option: use the paper as reference 0.4X0, X0=3.5mm
        // It was weight_masked += hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(0.4*3.5);
        // New version ?

        for (int j = 0; j < hit_energy->size(); j++) {
            if( masked && hit_isMasked->at(j) == 1 ) continue;
                nhit += 1;
                sume += hit_energy->at(j);
                weight += hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(3.5);
        }
    }
    return;
}

void GetPIDParametersTreeProcessor::hits_layer(float hlv[NUMBER_OF_LAYER], vector<float> * hit_energy, vector<int> *hit_slab, TVectorD W_thicknesses, vector<int> *hit_isMasked, bool masked=false, bool normalized=false, string count_type="nhit"){
  float hit_count[NUMBER_OF_LAYER];
  float sume[NUMBER_OF_LAYER];
  float sume_w[NUMBER_OF_LAYER];
  
  float sume_total = 1.;
  float sume_w_total = 1.;

  float weight = 1.;

  // Initial values
  for (int ilayer=0; ilayer < NUMBER_OF_LAYER; ilayer++){
    hit_count[ilayer] = 0. ;
    sume[ilayer] = 0. ;
    sume_w[ilayer] = 0. ;
  }
  
  if(hit_energy->size() > 0){
    for (int j = 0; j < hit_energy->size(); j++) {
      if( masked && hit_isMasked->at(j) == 1 ) continue;
      sume_total += hit_energy->at(j);
      sume_w_total += hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(3.5);
      //cout<<"|DEBUG(hit_layer)| sume total: "<<sume_total<<", sume_w_total: "<<sume_w_total<<endl;
      //cout<<"layer and thickness"<<hit_slab->at(j)<<" "<<W_thicknesses[hit_slab->at(j)]<<endl;
    }
  }
  
  for (int ilayer=0; ilayer < NUMBER_OF_LAYER; ilayer++) {
    if(hit_slab->size() > 0){
      for( int j = 0; j < hit_energy->size(); j++ ) {
	if( masked && hit_isMasked->at(j) == 1 ) continue;
	if( hit_slab->at(j) == ilayer ) {
	  hit_count[ilayer] += 1;
	  sume[ilayer] += hit_energy->at(j);
	  sume_w[ilayer] += hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(3.5);
	}
      }
      if(normalized == true) {
	if(count_type == "nhit") weight = 1./hit_energy->size();
	if(count_type == "sume") weight = 1./sume_total;
	if(count_type == "weight") weight = 1./sume_w_total;
      }  
      
      if(count_type == "nhit") hlv[ilayer] = hit_count[ilayer]*weight;
      if(count_type == "sume") hlv[ilayer] = sume[ilayer]*weight;
      if(count_type == "weight") hlv[ilayer] = sume_w[ilayer]*weight;
    }
    else hlv[ilayer] = 0.;
  
  }
  //cout<<"| DEBUG(hits_layer) | "<<"Mode: "<<count_type<<", normalized: "<<normalized<<endl;
  //cout<<"| DEBUG(hits_layer) | "<<"Hit energy size: "<<hit_energy->size()<<". Hit slab size: "<<hit_slab->size()<<endl;
  //cout<<"| DEBUG(hits_layer) | "<<"hlv: "<<hlv[0]<<" "<<hlv[1]<<" "<<hlv[2]<<" "<<hlv[3]<<" "<<hlv[4]<<" "<<hlv[5]<<" "<<hlv[6]<<" "<<hlv[7]<<" "<<hlv[8]<<" "<<hlv[9]<<" "<<hlv[10]<<" "<<hlv[11]<<" "<<hlv[12]<<" "<<hlv[13]<<" "<<hlv[14]<<endl;
    
  
  return ;
}

bool GetPIDParametersTreeProcessor::is_Shower(float entries, float array[NUMBER_OF_LAYER]) {
  float threshold = 3.;
  float shower_maxvalue = 0.;
  bool isShower = false;


  if(entries > 0){
    for(int ilayer=0; ilayer<NUMBER_OF_LAYER; ilayer++){
      float thislayer = array[ilayer];
      if(thislayer > shower_maxvalue){
        shower_maxvalue = thislayer;
      }
    }
    for(int ilayer=0; ilayer<NUMBER_OF_LAYER-2; ilayer++){
      float thislayer = array[ilayer];
      float nextlayer = array[ilayer+1];
      float nextnextlayer = array[ilayer+2];
      if((thislayer > threshold) and (nextlayer > thislayer) and (nextnextlayer > nextlayer) and (shower_maxvalue > 5.)){
        isShower = true;
        break;
      }
    }
  }
  return isShower;

}

void GetPIDParametersTreeProcessor::shower_variables(float entries, float array[NUMBER_OF_LAYER], float array_n[NUMBER_OF_LAYER], float &shower_maxvalue, float &shower_maxvalue_n, int &ilayermax,int &ilayerstart, int &ilayerstart_10, int &ilayerend, int &ilayerend_10, string count_type = "nhit", bool is_shower=false) {

  float percentage = 0.1;
  float threshold = 3.;
 
  if((entries > 0) and (is_shower == true)){
    for(int ilayer=0; ilayer<NUMBER_OF_LAYER; ilayer++){
      float thislayer = array[ilayer];
      float thislayer_n = array_n[ilayer];
      if(thislayer > shower_maxvalue){
        shower_maxvalue = thislayer;
        shower_maxvalue_n = thislayer_n;
        ilayermax = ilayer;
      }
    }
    for(int ilayer=NUMBER_OF_LAYER-1; ilayer>ilayermax; ilayer--){
      float thislayer = array[ilayer];
      if(thislayer > threshold){
        ilayerend = ilayer;
        break;
      }
    }
    for(int ilayer=NUMBER_OF_LAYER-1; ilayer>ilayermax; ilayer--){
      float thislayer = array[ilayer];
      if((thislayer > threshold) and (thislayer > percentage*shower_maxvalue)) {
        ilayerend_10 = ilayer;
        break;
      }
    }
 
    for(int ilayer=0; ilayer<ilayermax; ilayer++){
      float thislayer = array[ilayer];
       if((thislayer > threshold)) {
        ilayerstart = ilayer;
        break;
      }
    }
    for(int ilayer=0; ilayer<ilayermax; ilayer++){
      float thislayer = array[ilayer];
       if((thislayer > threshold) && (thislayer > percentage*shower_maxvalue)) {
        ilayerstart_10 = ilayer;
        break;
      }
    }

  }
  return;
}

float GetPIDParametersTreeProcessor::MIP_Likeness(float nhits_layer[NUMBER_OF_LAYER]) {
  float score = 0.;
  for(int i=0; i<NUMBER_OF_LAYER; i++) {
    if(nhits_layer[i] == 0) score -= 1.;
    if(nhits_layer[i] > 0) score += 1./(nhits_layer[i]);
  }
  score = (score/NUMBER_OF_LAYER + 1.)/2.;
  return score;
}

void GetPIDParametersTreeProcessor::is_interaction(float &ecal_int, int nhit_e) {
  ecal_int = 0.;
  if(nhit_e > 0) ecal_int = 1.;
  return;
}



float GetPIDParametersTreeProcessor::hits_max_distance(vector<int> *hit_slab, vector<float> * hit_x, vector<float> * hit_y, vector<int> * hit_isMasked, bool masked=false) {
  float max_distance = 0.;
  if(hit_slab->size() > 1){
    for (int i = 0; i < hit_slab->size(); i++){
      if (masked && hit_isMasked->at(i) == 1) continue;
      int layer1 = hit_slab->at(i);
      float hit1_x = hit_x->at(i);
      float hit1_y = hit_y->at(i);
      for (int j = i+1; j < hit_slab->size()-1; j++){
	if (masked && hit_isMasked->at(j) == 1) continue;
        int layer2 = hit_slab->at(j);
	if(layer2 != layer1) continue;
        float hit2_x = hit_x->at(j);
        float hit2_y = hit_y->at(j);
        float distance = pow(pow(hit2_x-hit1_x,2)+pow(hit2_y-hit1_y,2),0.5);
        if(distance>max_distance) max_distance = distance;
      }
    }
  }
  return max_distance;
}

/*struct GetPIDParametersTreeProcessor::hitpair
{
  float hit_rs;
  float hit_es;
};*/

bool GetPIDParametersTreeProcessor::CompareHitsR(const hitpair hit1, const hitpair hit2) { return hit1.hit_rs < hit2.hit_rs; }

float GetPIDParametersTreeProcessor::moliere(vector<float> * hit_energy, vector<int> *hit_slab, TVectorD W_thicknesses, vector<float> * hit_x, vector<float> * hit_y, vector<float> * hit_z,vector<int> * hit_isMasked, bool masked=false, float containment = 0.90, bool is_shower=false) {
  float mol_rad = 0.;
  float weighte = 0.;
  float wx = 0.; float wy = 0.; float wz = 0.;
  float r = 0.;

  if((hit_energy->size() > 0) and (is_shower==true)){
    vector<hitpair> hits_vector;

    for (int j = 0; j < hit_energy->size(); j++){
      if (masked && hit_isMasked->at(j) == 1) continue;
      weighte += hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(3.5);
      wx += hit_x->at(j) * hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(3.5);
      wy += hit_y->at(j) * hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(3.5);
    }

    float bary_x = wx / weighte;
    float bary_y = wy / weighte;
    if( (weighte == 0) or (weighte < 0) ){
      bary_x = 0.;
      bary_y = 0.;
    }

    for (int j = 0; j < hit_energy->size(); j++){
      if (masked && hit_isMasked->at(j) == 1) continue;
      r = pow(pow((hit_x->at(j) - bary_x) , 2) + pow((hit_y->at(j) - bary_y), 2), 0.5);
      hits_vector.push_back({static_cast<float>(r),static_cast<float>(hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(3.5))});
    }

    float mol_e = 0.;
    int mol_i=0;
    if(hits_vector.size() > 0){
      std::sort(hits_vector.begin(), hits_vector.end(), CompareHitsR);
      for (int j = 0; j < hits_vector.size(); j++) {
        mol_e += hits_vector.at(j).hit_es;
        if (mol_e >= containment * weighte){
          mol_i=j;
          break;
        }
      }

      if(mol_i<0) mol_i=0;
      mol_rad=hits_vector.at(mol_i).hit_rs;
    }
  }
  return mol_rad;
}

void GetPIDParametersTreeProcessor::radius_layer(float mol_per_layer[NUMBER_OF_LAYER], vector<float> * hit_energy, vector<int> *hit_slab, TVectorD W_thicknesses,vector<float> * hit_x, vector<float> * hit_y, vector<float> * hit_z,vector<int> * hit_isMasked, bool masked=false, float containment = 0.90, bool is_shower=false)  {

  if((hit_energy->size() > 0) and (is_shower==true)){

    for(int ilayer = 0; ilayer < NUMBER_OF_LAYER; ilayer++){
      float mol_rad = 0.;
      float weighte = 0.;
      float wx = 0.; float wy = 0.;
      float r = 0.;
      int nhit_layer = 0;
      struct hitcontent
      {
        float hit_rs;
        float hit_es;
      };
      vector<hitcontent> hits_vector;

      for (int j = 0; j < hit_energy->size(); j++){
        if (masked && hit_isMasked->at(j) == 1) continue;
        if(hit_slab->at(j) == ilayer) {
          weighte += hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(3.5);
          wx += hit_x->at(j) * hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(3.5);
          wy += hit_y->at(j) * hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(3.5);
          nhit_layer += 1;
        }
      }

      float bary_x = wx / weighte;
      float bary_y = wy / weighte;
      if( (weighte == 0) or (weighte < 0) ){
        bary_x = 0.;
        bary_y = 0.;
      }

      for (int j = 0; j < hit_energy->size(); j++){
        if (masked && hit_isMasked->at(j) == 1) continue;
        if (hit_slab->at(j) == ilayer) {
          r = pow(pow((hit_x->at(j) - bary_x) , 2) + pow((hit_y->at(j) - bary_y), 2), 0.5);
          hits_vector.push_back({static_cast<float>(r),static_cast<float>(hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(3.5))});
        }
      }
      float mol_e = 0.;
      int mol_i=0;
      if(hits_vector.size() > 0){

	std::sort(hits_vector.begin(), hits_vector.end(),
                  [](const auto& i, const auto& j) { return i.hit_rs < j.hit_rs; } );

        for (int j = 0; j < hits_vector.size(); j++) {
          mol_e += hits_vector.at(j).hit_es;
          if (mol_e >= containment * weighte){
            mol_i=j;
            break;
          }
        }

        if(mol_i<0) mol_i=0.;
        mol_per_layer[ilayer] = hits_vector.at(mol_i).hit_rs;
      }
      if(nhit_layer < 3) mol_per_layer[ilayer] = 0.;
      if(mol_per_layer[ilayer] < 1) mol_per_layer[ilayer] = 0.;
      hits_vector.clear();
    }
  }
  /*
  for(int iout = 0; iout < NUMBER_OF_LAYER; iout++) cout<<mol_per_layer[iout]<<" ";
  cout<<endl;
  */
  return;
}

void GetPIDParametersTreeProcessor::barycenter(vector<float> * hit_energy, vector<int> *hit_slab, TVectorD W_thicknesses,vector<float> * hit_x, vector<float> * hit_y, vector<float> * hit_z, float bar_xyzr[4],vector<int> * hit_isMasked, bool masked=false, bool is_shower=false) {

  float sume = 0.;
  float wx = 0.; float wy = 0.; float wz = 0.;
  float bary_x = 0., bary_y = 0., bary_z = 0.;
  // Removing shower condition
  //if((hit_energy->size() > 0) and (is_shower == true)){
  if(hit_energy->size() > 0){
    for (int j = 0; j < hit_energy->size(); j++){
      if (masked && hit_isMasked->at(j) == 1) continue;
      //cout<<"DBG "<<j<<endl;                                                                                                                                                                            
      sume += hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(3.5);
      wx += hit_x->at(j) * hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(3.5);
      wy += hit_y->at(j) * hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(3.5);
      wz += hit_z->at(j) * hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(3.5);
    }

    bary_x = wx / sume;
    bary_y = wy / sume;
    bary_z = wz / sume;

    if( (sume == 0) or (sume < 0) ){
      bary_x = 0.;
      bary_y = 0.;
      bary_z = 0.;
    }
  }

  bar_xyzr[0] = bary_x;
  bar_xyzr[1] = bary_y;
  bar_xyzr[2] = bary_z;
  bar_xyzr[3] = pow(pow(bary_x,2)+pow(bary_y,2),0.5);

  return;
}

void GetPIDParametersTreeProcessor::bary_layer(float blv[NUMBER_OF_LAYER][3], vector<float> * hit_energy, vector<int> *hit_slab,TVectorD W_thicknesses, vector<float> * hit_x, vector<float> * hit_y,vector<float> * hit_z, vector<int> * hit_isMasked, bool masked=false, bool is_shower=false) {

  float sume_w[NUMBER_OF_LAYER];
  float wx[NUMBER_OF_LAYER]; float wy[NUMBER_OF_LAYER];
  float bary_x[NUMBER_OF_LAYER]; float bary_y[NUMBER_OF_LAYER];

  for (int ilayer=0; ilayer < NUMBER_OF_LAYER; ilayer++){
    sume_w[ilayer] = 0. ;
    wx[ilayer] = 0.;
    wy[ilayer] = 0.;
    bary_x[ilayer] = 0.;
    bary_y[ilayer] = 0.;
    blv[ilayer][0] = 0.;
    blv[ilayer][1] = 0.;
    blv[ilayer][2] = 0.;
  }
  // Removing the shower condition
  // if((hit_energy->size() > 0) and (is_shower == true)){
  if(hit_energy->size() > 0){
    for (int ilayer=0; ilayer < NUMBER_OF_LAYER; ilayer++) {
      if(hit_slab->size() > 0){
        for (int j = 0; j < hit_energy->size(); j++) {
          if (masked && hit_isMasked->at(j) == 1) continue;
          if (hit_slab->at(j) == 0) {
            sume_w[ilayer] += hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(3.5);
            wx[ilayer] += hit_x->at(j) * hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(3.5);
            wy[ilayer] += hit_y->at(j) * hit_energy->at(j) * W_thicknesses[hit_slab->at(j)]/(3.5);
          }
        }
        bary_x[ilayer] = wx[ilayer] / sume_w[ilayer];
        bary_y[ilayer] = wy[ilayer] / sume_w[ilayer];
      }
      if( (sume_w[ilayer] == 0) or (sume_w[ilayer] < 0) ){
        bary_x[ilayer] = 0.;
        bary_y[ilayer] = 0.;
      }
      blv[ilayer][0] = bary_x[ilayer];
      blv[ilayer][1] = bary_y[ilayer];
      blv[ilayer][2] = pow(pow(bary_x[ilayer],2) + pow(bary_y[ilayer],2),0.5);
    }
  }

  return;
}

void GetPIDParametersTreeProcessor::graph_setup_add(TGraph *g, string title, Color_t color){
    g->SetTitle(title.c_str());
    g->SetLineColor(color);
    g->SetLineWidth(3);
    return;
}



void GetPIDParametersTreeProcessor::init() {
	AIDAProcessor::tree(this);
    //Si: 650 650 650 650 500 500 500 500 500 500 320 320 320 320 320
  
    TVectorD W_thicknesses(NUMBER_OF_LAYER, W);
   
   
    TTree *outtree = new TTree("ntp","NTuples");

    // branches definitions
    // We will save both the branches and the histos
    // Values

    // Adding the branches
    outtree->Branch("ecal_interaction",&b_ecal_interaction,"b_ecal_interaction/F");
    outtree->Branch("nhit",&b_nhit,"b_nhit/F");
    outtree->Branch("sume",&b_sume,"b_sume/F");
    outtree->Branch("weighte",&b_weighte,"b_weighte/F");
    outtree->Branch("mol",&b_mol,"mol/F");
    outtree->Branch("MIP_Likeness",&b_MIP_Likeness,"MIP_Likeness/F");
    outtree->Branch("hits_max_distance",&b_hits_max_distance,"b_hits_max_distance/F");
    outtree->Branch("radius90_layer_0",&b_radius90_layer_0,"b_radius90_layer_0/F");
    outtree->Branch("radius90_layer_1",&b_radius90_layer_1,"b_radius90_layer_1/F");
    outtree->Branch("radius90_layer_2",&b_radius90_layer_2,"b_radius90_layer_2/F");
    outtree->Branch("radius90_layer_3",&b_radius90_layer_3,"b_radius90_layer_3/F");
    outtree->Branch("radius90_layer_4",&b_radius90_layer_4,"b_radius90_layer_4/F");
    outtree->Branch("radius90_layer_5",&b_radius90_layer_5,"b_radius90_layer_5/F");
    outtree->Branch("radius90_layer_6",&b_radius90_layer_6,"b_radius90_layer_6/F");
    outtree->Branch("radius90_layer_7",&b_radius90_layer_7,"b_radius90_layer_7/F");
    outtree->Branch("radius90_layer_8",&b_radius90_layer_8,"b_radius90_layer_8/F");
    outtree->Branch("radius90_layer_9",&b_radius90_layer_9,"b_radius90_layer_9/F");
    outtree->Branch("radius90_layer_10",&b_radius90_layer_10,"b_radius90_layer_10/F");
    outtree->Branch("radius90_layer_11",&b_radius90_layer_11,"b_radius90_layer_11/F");
    outtree->Branch("radius90_layer_12",&b_radius90_layer_12,"b_radius90_layer_12/F");
    outtree->Branch("radius90_layer_13",&b_radius90_layer_13,"b_radius90_layer_13/F");
    outtree->Branch("radius90_layer_14",&b_radius90_layer_14,"b_radius90_layer_14/F");
    outtree->Branch("bar_x",&b_bar_x,"b_bar_x/F");
    outtree->Branch("bar_y",&b_bar_y,"b_bar_y/F");
    outtree->Branch("bar_z",&b_bar_z,"b_bar_z/F");
    outtree->Branch("bar_r",&b_bar_r,"b_bar_r/F");
    outtree->Branch("bar_x_layer_0",&b_bar_x_layer_0,"b_bar_x_layer_0/F");
    outtree->Branch("bar_x_layer_1",&b_bar_x_layer_1,"b_bar_x_layer_1/F");
    outtree->Branch("bar_x_layer_2",&b_bar_x_layer_2,"b_bar_x_layer_2/F");
    outtree->Branch("bar_x_layer_3",&b_bar_x_layer_3,"b_bar_x_layer_3/F");
    outtree->Branch("bar_x_layer_4",&b_bar_x_layer_4,"b_bar_x_layer_4/F");
    outtree->Branch("bar_x_layer_5",&b_bar_x_layer_5,"b_bar_x_layer_5/F");
    outtree->Branch("bar_x_layer_6",&b_bar_x_layer_6,"b_bar_x_layer_6/F");
    outtree->Branch("bar_x_layer_7",&b_bar_x_layer_7,"b_bar_x_layer_7/F");
    outtree->Branch("bar_x_layer_8",&b_bar_x_layer_8,"b_bar_x_layer_8/F");
    outtree->Branch("bar_x_layer_9",&b_bar_x_layer_9,"b_bar_x_layer_9/F");
    outtree->Branch("bar_x_layer_10",&b_bar_x_layer_10,"b_bar_x_layer_10/F");
    outtree->Branch("bar_x_layer_11",&b_bar_x_layer_11,"b_bar_x_layer_11/F");
    outtree->Branch("bar_x_layer_12",&b_bar_x_layer_12,"b_bar_x_layer_12/F");
    outtree->Branch("bar_x_layer_13",&b_bar_x_layer_13,"b_bar_x_layer_13/F");
    outtree->Branch("bar_x_layer_14",&b_bar_x_layer_14,"b_bar_x_layer_14/F");
    outtree->Branch("bar_y_layer_0",&b_bar_y_layer_0,"b_bar_y_layer_0/F");
    outtree->Branch("bar_y_layer_1",&b_bar_y_layer_1,"b_bar_y_layer_1/F");
    outtree->Branch("bar_y_layer_2",&b_bar_y_layer_2,"b_bar_y_layer_2/F");
    outtree->Branch("bar_y_layer_3",&b_bar_y_layer_3,"b_bar_y_layer_3/F");
    outtree->Branch("bar_y_layer_4",&b_bar_y_layer_4,"b_bar_y_layer_4/F");
    outtree->Branch("bar_y_layer_5",&b_bar_y_layer_5,"b_bar_y_layer_5/F");
    outtree->Branch("bar_y_layer_6",&b_bar_y_layer_6,"b_bar_y_layer_6/F");
    outtree->Branch("bar_y_layer_7",&b_bar_y_layer_7,"b_bar_y_layer_7/F");
    outtree->Branch("bar_y_layer_8",&b_bar_y_layer_8,"b_bar_y_layer_8/F");
    outtree->Branch("bar_y_layer_9",&b_bar_y_layer_9,"b_bar_y_layer_9/F");
    outtree->Branch("bar_y_layer_10",&b_bar_y_layer_10,"b_bar_y_layer_10/F");
    outtree->Branch("bar_y_layer_11",&b_bar_y_layer_11,"b_bar_y_layer_11/F");
    outtree->Branch("bar_y_layer_12",&b_bar_y_layer_12,"b_bar_y_layer_12/F");
    outtree->Branch("bar_y_layer_13",&b_bar_y_layer_13,"b_bar_y_layer_13/F");
    outtree->Branch("bar_y_layer_14",&b_bar_y_layer_14,"b_bar_y_layer_14/F");
    outtree->Branch("bar_r_layer_0",&b_bar_r_layer_0,"b_bar_r_layer_0/F");
    outtree->Branch("bar_r_layer_1",&b_bar_r_layer_1,"b_bar_r_layer_1/F");
    outtree->Branch("bar_r_layer_2",&b_bar_r_layer_2,"b_bar_r_layer_2/F");
    outtree->Branch("bar_r_layer_3",&b_bar_r_layer_3,"b_bar_r_layer_3/F");
    outtree->Branch("bar_r_layer_4",&b_bar_r_layer_4,"b_bar_r_layer_4/F");
    outtree->Branch("bar_r_layer_5",&b_bar_r_layer_5,"b_bar_r_layer_5/F");
    outtree->Branch("bar_r_layer_6",&b_bar_r_layer_6,"b_bar_r_layer_6/F");
    outtree->Branch("bar_r_layer_7",&b_bar_r_layer_7,"b_bar_r_layer_7/F");
    outtree->Branch("bar_r_layer_8",&b_bar_r_layer_8,"b_bar_r_layer_8/F");
    outtree->Branch("bar_r_layer_9",&b_bar_r_layer_9,"b_bar_r_layer_9/F");
    outtree->Branch("bar_r_layer_10",&b_bar_r_layer_10,"b_bar_r_layer_10/F");
    outtree->Branch("bar_r_layer_11",&b_bar_r_layer_11,"b_bar_r_layer_11/F");
    outtree->Branch("bar_r_layer_12",&b_bar_r_layer_12,"b_bar_r_layer_12/F");
    outtree->Branch("bar_r_layer_13",&b_bar_r_layer_13,"b_bar_r_layer_13/F");
    outtree->Branch("bar_r_layer_14",&b_bar_r_layer_14,"b_bar_r_layer_14/F");
    outtree->Branch("shower_nhit_max_layer",&b_shower_nhit_max_layer,"b_shower_nhit_max_layer/F");
    outtree->Branch("shower_nhit_start_layer",&b_shower_nhit_start_layer,"b_shower_nhit_start_layer/F");
    outtree->Branch("shower_nhit_end_layer",&b_shower_nhit_end_layer,"b_shower_nhit_end_layer/F");
    outtree->Branch("shower_nhit_start_10_layer",&b_shower_nhit_start_10_layer,"b_shower_nhit_start_10_layer/F");
    outtree->Branch("shower_nhit_end_10_layer",&b_shower_nhit_end_10_layer,"b_shower_nhit_end_10_layer/F");
    outtree->Branch("shower_nhit_average",&b_shower_nhit_average,"b_shower_nhit_average/F");
    outtree->Branch("shower_nhit_max",&b_shower_nhit_max,"b_shower_nhit_max/F");
    outtree->Branch("shower_sume_max_layer",&b_shower_sume_max_layer,"b_shower_sume_max_layer/F");
    outtree->Branch("shower_sume_start_layer",&b_shower_sume_start_layer,"b_shower_sume_start_layer/F");
    outtree->Branch("shower_sume_end_layer",&b_shower_sume_end_layer,"b_shower_sume_end_layer/F");
    outtree->Branch("shower_sume_start_10_layer",&b_shower_sume_start_10_layer,"b_shower_sume_start_10_layer/F");
    outtree->Branch("shower_sume_end_10_layer",&b_shower_sume_end_10_layer,"b_shower_sume_end_10_layer/F");
    outtree->Branch("shower_sume_average",&b_shower_sume_average,"b_shower_sume_average/F");
    outtree->Branch("shower_sume_max",&b_shower_sume_max,"b_shower_sume_max/F");
    outtree->Branch("shower_weighte_max_layer",&b_shower_weighte_max_layer,"b_shower_weighte_max_layer/F");
    outtree->Branch("shower_weighte_start_layer",&b_shower_weighte_start_layer,"b_shower_weighte_start_layer/F");
    outtree->Branch("shower_weighte_end_layer",&b_shower_weighte_end_layer,"b_shower_weighte_end_layer/F");
    outtree->Branch("shower_weighte_start_10_layer",&b_shower_weighte_start_10_layer,"b_shower_weighte_start_10_layer/F");
    outtree->Branch("shower_weighte_end_10_layer",&b_shower_weighte_end_10_layer,"b_shower_weighte_end_10_layer/F");
    outtree->Branch("shower_weighte_average",&b_shower_weighte_average,"b_shower_weighte_average/F");
    outtree->Branch("shower_weighte_max",&b_shower_weighte_max,"b_shower_weighte_max/F");
    outtree->Branch("nhit_layer_0",&b_nhit_layer_0,"b_nhit_layer_0/F");
    outtree->Branch("nhit_layer_1",&b_nhit_layer_1,"b_nhit_layer_1/F");
    outtree->Branch("nhit_layer_2",&b_nhit_layer_2,"b_nhit_layer_2/F");
    outtree->Branch("nhit_layer_3",&b_nhit_layer_3,"b_nhit_layer_3/F");
    outtree->Branch("nhit_layer_4",&b_nhit_layer_4,"b_nhit_layer_4/F");
    outtree->Branch("nhit_layer_5",&b_nhit_layer_5,"b_nhit_layer_5/F");
    outtree->Branch("nhit_layer_6",&b_nhit_layer_6,"b_nhit_layer_6/F");
    outtree->Branch("nhit_layer_7",&b_nhit_layer_7,"b_nhit_layer_7/F");
    outtree->Branch("nhit_layer_8",&b_nhit_layer_8,"b_nhit_layer_8/F");
    outtree->Branch("nhit_layer_9",&b_nhit_layer_9,"b_nhit_layer_9/F");
    outtree->Branch("nhit_layer_10",&b_nhit_layer_10,"b_nhit_layer_10/F");
    outtree->Branch("nhit_layer_11",&b_nhit_layer_11,"b_nhit_layer_11/F");
    outtree->Branch("nhit_layer_12",&b_nhit_layer_12,"b_nhit_layer_12/F");
    outtree->Branch("nhit_layer_13",&b_nhit_layer_13,"b_nhit_layer_13/F");
    outtree->Branch("nhit_layer_14",&b_nhit_layer_14,"b_nhit_layer_14/F");
    outtree->Branch("nhit_layer_n_0",&b_nhit_layer_n_0,"b_nhit_layer_n_0/F");
    outtree->Branch("nhit_layer_n_1",&b_nhit_layer_n_1,"b_nhit_layer_n_1/F");
    outtree->Branch("nhit_layer_n_2",&b_nhit_layer_n_2,"b_nhit_layer_n_2/F");
    outtree->Branch("nhit_layer_n_3",&b_nhit_layer_n_3,"b_nhit_layer_n_3/F");
    outtree->Branch("nhit_layer_n_4",&b_nhit_layer_n_4,"b_nhit_layer_n_4/F");
    outtree->Branch("nhit_layer_n_5",&b_nhit_layer_n_5,"b_nhit_layer_n_5/F");
    outtree->Branch("nhit_layer_n_6",&b_nhit_layer_n_6,"b_nhit_layer_n_6/F");
    outtree->Branch("nhit_layer_n_7",&b_nhit_layer_n_7,"b_nhit_layer_n_7/F");
    outtree->Branch("nhit_layer_n_8",&b_nhit_layer_n_8,"b_nhit_layer_n_8/F");
    outtree->Branch("nhit_layer_n_9",&b_nhit_layer_n_9,"b_nhit_layer_n_9/F");
    outtree->Branch("nhit_layer_n_10",&b_nhit_layer_n_10,"b_nhit_layer_n_10/F");
    outtree->Branch("nhit_layer_n_11",&b_nhit_layer_n_11,"b_nhit_layer_n_11/F");
    outtree->Branch("nhit_layer_n_12",&b_nhit_layer_n_12,"b_nhit_layer_n_12/F");
    outtree->Branch("nhit_layer_n_13",&b_nhit_layer_n_13,"b_nhit_layer_n_13/F");
    outtree->Branch("nhit_layer_n_14",&b_nhit_layer_n_14,"b_nhit_layer_n_14/F");
    outtree->Branch("sume_layer_0",&b_sume_layer_0,"b_sume_layer_0/F");
    outtree->Branch("sume_layer_1",&b_sume_layer_1,"b_sume_layer_1/F");
    outtree->Branch("sume_layer_2",&b_sume_layer_2,"b_sume_layer_2/F");
    outtree->Branch("sume_layer_3",&b_sume_layer_3,"b_sume_layer_3/F");
    outtree->Branch("sume_layer_4",&b_sume_layer_4,"b_sume_layer_4/F");
    outtree->Branch("sume_layer_5",&b_sume_layer_5,"b_sume_layer_5/F");
    outtree->Branch("sume_layer_6",&b_sume_layer_6,"b_sume_layer_6/F");
    outtree->Branch("sume_layer_7",&b_sume_layer_7,"b_sume_layer_7/F");
    outtree->Branch("sume_layer_8",&b_sume_layer_8,"b_sume_layer_8/F");
    outtree->Branch("sume_layer_9",&b_sume_layer_9,"b_sume_layer_9/F");
    outtree->Branch("sume_layer_10",&b_sume_layer_10,"b_sume_layer_10/F");
    outtree->Branch("sume_layer_11",&b_sume_layer_11,"b_sume_layer_11/F");
    outtree->Branch("sume_layer_12",&b_sume_layer_12,"b_sume_layer_12/F");
    outtree->Branch("sume_layer_13",&b_sume_layer_13,"b_sume_layer_13/F");
    outtree->Branch("sume_layer_14",&b_sume_layer_14,"b_sume_layer_14/F");
    outtree->Branch("sume_layer_n_0",&b_sume_layer_n_0,"b_sume_layer_n_0/F");
    outtree->Branch("sume_layer_n_1",&b_sume_layer_n_1,"b_sume_layer_n_1/F");
    outtree->Branch("sume_layer_n_2",&b_sume_layer_n_2,"b_sume_layer_n_2/F");
    outtree->Branch("sume_layer_n_3",&b_sume_layer_n_3,"b_sume_layer_n_3/F");
    outtree->Branch("sume_layer_n_4",&b_sume_layer_n_4,"b_sume_layer_n_4/F");
    outtree->Branch("sume_layer_n_5",&b_sume_layer_n_5,"b_sume_layer_n_5/F");
    outtree->Branch("sume_layer_n_6",&b_sume_layer_n_6,"b_sume_layer_n_6/F");
    outtree->Branch("sume_layer_n_7",&b_sume_layer_n_7,"b_sume_layer_n_7/F");
    outtree->Branch("sume_layer_n_8",&b_sume_layer_n_8,"b_sume_layer_n_8/F");
    outtree->Branch("sume_layer_n_9",&b_sume_layer_n_9,"b_sume_layer_n_9/F");
    outtree->Branch("sume_layer_n_10",&b_sume_layer_n_10,"b_sume_layer_n_10/F");
    outtree->Branch("sume_layer_n_11",&b_sume_layer_n_11,"b_sume_layer_n_11/F");
    outtree->Branch("sume_layer_n_12",&b_sume_layer_n_12,"b_sume_layer_n_12/F");
    outtree->Branch("sume_layer_n_13",&b_sume_layer_n_13,"b_sume_layer_n_13/F");
    outtree->Branch("sume_layer_n_14",&b_sume_layer_n_14,"b_sume_layer_n_14/F");
    outtree->Branch("weighte_layer_0",&b_weighte_layer_0,"b_weighte_layer_0/F");
    outtree->Branch("weighte_layer_1",&b_weighte_layer_1,"b_weighte_layer_1/F");
    outtree->Branch("weighte_layer_2",&b_weighte_layer_2,"b_weighte_layer_2/F");
    outtree->Branch("weighte_layer_3",&b_weighte_layer_3,"b_weighte_layer_3/F");
    outtree->Branch("weighte_layer_4",&b_weighte_layer_4,"b_weighte_layer_4/F");
    outtree->Branch("weighte_layer_5",&b_weighte_layer_5,"b_weighte_layer_5/F");
    outtree->Branch("weighte_layer_6",&b_weighte_layer_6,"b_weighte_layer_6/F");
    outtree->Branch("weighte_layer_7",&b_weighte_layer_7,"b_weighte_layer_7/F");
    outtree->Branch("weighte_layer_8",&b_weighte_layer_8,"b_weighte_layer_8/F");
    outtree->Branch("weighte_layer_9",&b_weighte_layer_9,"b_weighte_layer_9/F");
    outtree->Branch("weighte_layer_10",&b_weighte_layer_10,"b_weighte_layer_10/F");
    outtree->Branch("weighte_layer_11",&b_weighte_layer_11,"b_weighte_layer_11/F");
    outtree->Branch("weighte_layer_12",&b_weighte_layer_12,"b_weighte_layer_12/F");
    outtree->Branch("weighte_layer_13",&b_weighte_layer_13,"b_weighte_layer_13/F");
    outtree->Branch("weighte_layer_14",&b_weighte_layer_14,"b_weighte_layer_14/F");
    outtree->Branch("weighte_layer_n_0",&b_weighte_layer_n_0,"b_weighte_layer_n_0/F");
    outtree->Branch("weighte_layer_n_1",&b_weighte_layer_n_1,"b_weighte_layer_n_1/F");
    outtree->Branch("weighte_layer_n_2",&b_weighte_layer_n_2,"b_weighte_layer_n_2/F");
    outtree->Branch("weighte_layer_n_3",&b_weighte_layer_n_3,"b_weighte_layer_n_3/F");
    outtree->Branch("weighte_layer_n_4",&b_weighte_layer_n_4,"b_weighte_layer_n_4/F");
    outtree->Branch("weighte_layer_n_5",&b_weighte_layer_n_5,"b_weighte_layer_n_5/F");
    outtree->Branch("weighte_layer_n_6",&b_weighte_layer_n_6,"b_weighte_layer_n_6/F");
    outtree->Branch("weighte_layer_n_7",&b_weighte_layer_n_7,"b_weighte_layer_n_7/F");
    outtree->Branch("weighte_layer_n_8",&b_weighte_layer_n_8,"b_weighte_layer_n_8/F");
    outtree->Branch("weighte_layer_n_9",&b_weighte_layer_n_9,"b_weighte_layer_n_9/F");
    outtree->Branch("weighte_layer_n_10",&b_weighte_layer_n_10,"b_weighte_layer_n_10/F");
    outtree->Branch("weighte_layer_n_11",&b_weighte_layer_n_11,"b_weighte_layer_n_11/F");
    outtree->Branch("weighte_layer_n_12",&b_weighte_layer_n_12,"b_weighte_layer_n_12/F");
    outtree->Branch("weighte_layer_n_13",&b_weighte_layer_n_13,"b_weighte_layer_n_13/F");
    outtree->Branch("weighte_layer_n_14",&b_weighte_layer_n_14,"b_weighte_layer_n_14/F");
    outtree->Branch("sume_layer_0",&b_sume_layer_0,"b_sume_layer_0/F");
    outtree->Branch("sume_layer_1",&b_sume_layer_1,"b_sume_layer_1/F");
    outtree->Branch("sume_layer_2",&b_sume_layer_2,"b_sume_layer_2/F");
    outtree->Branch("sume_layer_3",&b_sume_layer_3,"b_sume_layer_3/F");
    outtree->Branch("sume_layer_4",&b_sume_layer_4,"b_sume_layer_4/F");
    outtree->Branch("sume_layer_5",&b_sume_layer_5,"b_sume_layer_5/F");
    outtree->Branch("sume_layer_6",&b_sume_layer_6,"b_sume_layer_6/F");
    outtree->Branch("sume_layer_7",&b_sume_layer_7,"b_sume_layer_7/F");
    outtree->Branch("sume_layer_8",&b_sume_layer_8,"b_sume_layer_8/F");
    outtree->Branch("sume_layer_9",&b_sume_layer_9,"b_sume_layer_9/F");
    outtree->Branch("sume_layer_10",&b_sume_layer_10,"b_sume_layer_10/F");
    outtree->Branch("sume_layer_11",&b_sume_layer_11,"b_sume_layer_11/F");
    outtree->Branch("sume_layer_12",&b_sume_layer_12,"b_sume_layer_12/F");
    outtree->Branch("sume_layer_13",&b_sume_layer_13,"b_sume_layer_13/F");
    outtree->Branch("sume_layer_14",&b_sume_layer_14,"b_sume_layer_14/F");
    outtree->Branch("sume_layer_n_0",&b_sume_layer_n_0,"b_sume_layer_n_0/F");
    outtree->Branch("sume_layer_n_1",&b_sume_layer_n_1,"b_sume_layer_n_1/F");
    outtree->Branch("sume_layer_n_2",&b_sume_layer_n_2,"b_sume_layer_n_2/F");
    outtree->Branch("sume_layer_n_3",&b_sume_layer_n_3,"b_sume_layer_n_3/F");
    outtree->Branch("sume_layer_n_4",&b_sume_layer_n_4,"b_sume_layer_n_4/F");
    outtree->Branch("sume_layer_n_5",&b_sume_layer_n_5,"b_sume_layer_n_5/F");
    outtree->Branch("sume_layer_n_6",&b_sume_layer_n_6,"b_sume_layer_n_6/F");
    outtree->Branch("sume_layer_n_7",&b_sume_layer_n_7,"b_sume_layer_n_7/F");
    outtree->Branch("sume_layer_n_8",&b_sume_layer_n_8,"b_sume_layer_n_8/F");
    outtree->Branch("sume_layer_n_9",&b_sume_layer_n_9,"b_sume_layer_n_9/F");
    outtree->Branch("sume_layer_n_10",&b_sume_layer_n_10,"b_sume_layer_n_10/F");
    outtree->Branch("sume_layer_n_11",&b_sume_layer_n_11,"b_sume_layer_n_11/F");
    outtree->Branch("sume_layer_n_12",&b_sume_layer_n_12,"b_sume_layer_n_12/F");
    outtree->Branch("sume_layer_n_13",&b_sume_layer_n_13,"b_sume_layer_n_13/F");
    outtree->Branch("sume_layer_n_14",&b_sume_layer_n_14,"b_sume_layer_n_14/F");
    // Resolution histos
    string part_string = "e";//particle;

	printParameters();
    if (_MCColName!="Not configured in xml file") {_flagMcCol = true;}
    if (_ECALColName!="Not configured in xml file") {_flagEcalCol = true;}
    if (_pECALColName!="Not configured in xml file") {_flagPixelEcalCol = true;}
    if (_dECALColName!="Not configured in xml file") {_flagDigitEcalCol = true;}
    // MC particle
    if (_flagMcCol) {
    }

    // Monolithic calorimeter hits
    if (_flagEcalCol) {
        for (int i = 0; i < NUMBER_OF_LAYER; i++) {
        }
    }

    // Pixelised calorimeter hits
    if (_flagPixelEcalCol) {
        for (int i = 0; i < NUMBER_OF_LAYER; i++) {
        }
    }

    // Digitised calorimeter hits
    if (_flagDigitEcalCol) {
        // The histogram will instead be declared and filled at the ending stage
        // Bin ranges will be changed at the final stage
        for (int i = 0; i < NUMBER_OF_LAYER; i++) {
        }
    }
}


void GetPIDParametersTreeProcessor::ShowMCInfo(EVENT::LCCollection *myCollection) {
    int number = myCollection->getNumberOfElements();
    
    if (runEnergy==-1) {
        for (int i = 0; i < number; i++) {//Loop through the MC Particle collection for one event
            MCParticle *particle = dynamic_cast<MCParticle *>(myCollection->getElementAt(i));
            vector<MCParticle *> daughters = particle->getDaughters();
            runEnergy = particle->getEnergy();
            
            streamlog_out(DEBUG) << "\n MCCollection, particle:" << i;
            streamlog_out(DEBUG) << " pdg = " << particle->getPDG() <<",";
            streamlog_out(DEBUG) << " status = " << particle->getGeneratorStatus() <<",";
            streamlog_out(DEBUG) << " N_daughters = " << daughters.size() <<",";
            streamlog_out(DEBUG) << " E = " << particle->getEnergy() <<" GeV,";
            streamlog_out(DEBUG) << " px = " << particle->getMomentum()[0] <<" GeV,";
            streamlog_out(DEBUG) << " py = " << particle->getMomentum()[1] <<" GeV,";
            streamlog_out(DEBUG) << " pz = " << particle->getMomentum()[2] <<" GeV,";
            streamlog_out(DEBUG) << " m = " << particle->getMass() <<" GeV,";
            streamlog_out(DEBUG) << " charge = " << particle->getCharge() <<".";
        }
        streamlog_out(DEBUG) << std::endl;
    }
}

void GetPIDParametersTreeProcessor::ShowECALInfo(EVENT::LCCollection *myCollection) {
    int number = myCollection->getNumberOfElements();
    streamlog_out(DEBUG) << "TOTAL NUMBER OF HITS: " << number <<endl;
    CellIDDecoder<EVENT::SimCalorimeterHit> cd(myCollection);
    
	double totalEnergy = 0;

    for (int i = 0; i < number; i++) {
        SimCalorimeterHit *ecalhit = dynamic_cast<SimCalorimeterHit *>(myCollection->getElementAt(i));

        int xyz_x = cd(ecalhit)["x"];
        int xyz_y = cd(ecalhit)["y"];
        int xyz_z = cd(ecalhit)["layer"];
        float hit_energy = ecalhit->getEnergy();

        streamlog_out(DEBUG) << "\n SimCalorimeterHit, :" << i;
        streamlog_out(DEBUG) << " cellID-encoded=" << ecalhit->getCellID0();
        streamlog_out(DEBUG) << " x = " << xyz_x <<" mm,";
        streamlog_out(DEBUG) << " y = " << xyz_y <<" mm,";
        streamlog_out(DEBUG) << " z = " << xyz_z <<" layer,";
        streamlog_out(DEBUG) << " energy = " << hit_energy <<"GeV.\n";
        totalEnergy += hit_energy;
    }
    // return totalEnergyLayerSi;
	for (int i = 0; i < NUMBER_OF_LAYER; i++) {
	}
}//By this point all histograms are filled for one event, this is repeated for all the events in the collection

void GetPIDParametersTreeProcessor::ShowPixelECALInfo(EVENT::LCCollection *myCollection) {
    
    TVectorD W_thicknesses(NUMBER_OF_LAYER, W);
    int number = myCollection->getNumberOfElements();
    streamlog_out(DEBUG) << "TOTAL NUMBER OF HITS: " << number <<endl;
    CellIDDecoder<EVENT::SimCalorimeterHit> cd(myCollection);
    
	double totalEnergy = 0;
    int totalHits = 0;
    
    vector<float> *hit_energy = 0;
    vector<float> *hit_x = 0;
    vector<float> *hit_y = 0;
    vector<float> *hit_z = 0;
    vector<int> *hit_isMasked = 0;
    vector<int> *hit_slab = 0;
    
    for (int i = 0; i < number; i++) {
        SimCalorimeterHit *ecalhit = dynamic_cast<SimCalorimeterHit *>(myCollection->getElementAt(i));
        
        int IJK_I = cd(ecalhit)["I"];
        int IJK_J = cd(ecalhit)["J"];
        int IJK_K = cd(ecalhit)["K"];
        float hit_energyf = ecalhit->getEnergy();

        streamlog_out(DEBUG) << "\n SimCalorimeterHit, :" << i;
        streamlog_out(DEBUG) << " cellID-encoded=" << ecalhit->getCellID0();
        streamlog_out(DEBUG) << " I = " << IJK_I <<" mm,";
        streamlog_out(DEBUG) << " J = " << IJK_J <<" mm,";
        streamlog_out(DEBUG) << " K = " << IJK_K <<" layer,";
        streamlog_out(DEBUG) << " energy = " << hit_energy <<" GeV.\n";
		//totalEnergy += hit_energy;
    //    totalHits++;
        
        hit_energy->push_back(hit_energyf);
        hit_x->push_back(IJK_I);
        hit_y->push_back(IJK_I);
        hit_z->push_back(IJK_I);
    }
        
        // Resolution
        int nhit = 0;     
        float sume = 0;   
        float weighte = 0; 
        float mol_value = 0.;
        float bar_xyzr[4]; // 4 for (x,y,z,r)
        
        float nhit_layer_array[NUMBER_OF_LAYER];
        float nhit_layer_n_array[NUMBER_OF_LAYER];

        float sume_layer_array[NUMBER_OF_LAYER];
            float sume_layer_n_array[NUMBER_OF_LAYER];

        float weighte_layer_array[NUMBER_OF_LAYER];
            float weighte_layer_n_array[NUMBER_OF_LAYER];

        float bar_layer_array[NUMBER_OF_LAYER][3]; // 3 for (x,y,r)
        float radius90_layer_array[NUMBER_OF_LAYER];

        for(int i = 0; i<NUMBER_OF_LAYER;i++) {
                nhit_layer_array[i] = 0.;
                nhit_layer_n_array[i] = 0.;
                sume_layer_array[i] = 0.;
                sume_layer_n_array[i] = 0.;
                weighte_layer_array[i] = 0.;
                weighte_layer_n_array[i] = 0.;
                bar_layer_array[i][0] = 0.;
                bar_layer_array[i][1] = 0.;
                bar_layer_array[i][2] = 0.;
                radius90_layer_array[i] = 0.;
            }
                
        
        get_res(nhit,	sume, weighte, hit_energy, hit_slab, W_thicknesses, hit_isMasked, masked);
        
        // Fill barycenter
        barycenter(hit_energy, hit_slab, W_thicknesses, hit_x, hit_y, hit_z, bar_xyzr, hit_isMasked, masked);
        
        // Fill shower profile 
        // Nhit
        hits_layer(nhit_layer_array, hit_energy, hit_slab, W_thicknesses, hit_isMasked, masked, false, "nhit");
        hits_layer(nhit_layer_n_array, hit_energy, hit_slab, W_thicknesses, hit_isMasked, masked, true, "nhit");
        
        // Fill MIP-Likeness
        float MIP_Likeness_value = MIP_Likeness(nhit_layer_array);

        // Sume
            hits_layer(sume_layer_array, hit_energy, hit_slab, W_thicknesses, hit_isMasked, masked, false, "sume");
            hits_layer(sume_layer_n_array, hit_energy, hit_slab, W_thicknesses, hit_isMasked, masked, true, "sume");

        // Weighted energy
            hits_layer(weighte_layer_array, hit_energy, hit_slab, W_thicknesses, hit_isMasked, masked, false, "weight");
            hits_layer(weighte_layer_n_array, hit_energy, hit_slab, W_thicknesses, hit_isMasked, masked, true, "weight");

        bool shower_bool = is_Shower(sume, sume_layer_array);

        is_interaction(b_ecal_interaction, nhit);
        
        
        // shower nhit general parameters start
        float nhit_shower_maxvalue = 0.;
        float nhit_shower_maxvalue_n = 0.;
        int nhit_ilayermax = -1;
        int nhit_ilayerstart = -1;
        int nhit_ilayerend = -1;
        int nhit_ilayerstart_10 = -1;
            int nhit_ilayerend_10 = -1;
        
        float nhit_shower_averagevalue = nhit/NUMBER_OF_LAYER;

        shower_variables(nhit, nhit_layer_array, nhit_layer_n_array, nhit_shower_maxvalue, nhit_shower_maxvalue_n, nhit_ilayermax,nhit_ilayerstart, nhit_ilayerstart_10, nhit_ilayerend, nhit_ilayerend_10, "nhit", shower_bool);
        //cout<<"nhit shower variables: maxvalue, maxvalue_n,ilayermax, ilayerstart, ilayerstart(10%), ilayerend, ilayerend(10%)"<<endl;
        //cout<<"nhit shower variables: "<<nhit<<" "<<nhit_shower_maxvalue<<" "<<nhit_shower_maxvalue_n<<" "<<nhit_ilayermax<<" "<<
        //  nhit_ilayerstart<<" "<<nhit_ilayerstart_10<<" "<<nhit_ilayerend<<" "<<nhit_ilayerend_10<<endl;

        //shower nhit general parameters finish
        
            // shower sume general parameters start
            float sume_shower_maxvalue = 0.;
            float sume_shower_maxvalue_n = 0.;
            int sume_ilayermax = -1;
            int sume_ilayerstart = -1;
            int sume_ilayerend = -1;
            int sume_ilayerstart_10 = -1;
            int sume_ilayerend_10 = -1;

            float sume_shower_averagevalue = sume/NUMBER_OF_LAYER;

            shower_variables(sume, sume_layer_array, sume_layer_n_array, sume_shower_maxvalue, sume_shower_maxvalue_n, sume_ilayermax,
                            sume_ilayerstart, sume_ilayerstart_10, sume_ilayerend, sume_ilayerend_10, "sume", shower_bool);
            
            //shower sume general parameters finish

        // shower weight general parameters start
            float weighte_shower_maxvalue = 0.;
        float weighte_shower_maxvalue_n = 0.;
            int weighte_ilayermax = -1;
            int weighte_ilayerstart = -1;
            int weighte_ilayerend = -1;
            int weighte_ilayerstart_10 = -1;
            int weighte_ilayerend_10 = -1;

        float weighte_shower_averagevalue = weighte/NUMBER_OF_LAYER;

        shower_variables(weighte, weighte_layer_array, weighte_layer_n_array, weighte_shower_maxvalue, weighte_shower_maxvalue_n, weighte_ilayermax,
                            weighte_ilayerstart, weighte_ilayerstart_10, weighte_ilayerend, weighte_ilayerend_10, "weight", shower_bool);
        //cout<<"weight shower variables: maxvalue, maxvalue_n,ilayermax, ilayerstart, ilayerstart(10%), ilayerend, ilayerend(10%)"<<endl;
        //cout<<"weight shower variables: "<<weight<<" "<<weighte_shower_maxvalue<<" "<<weighte_shower_maxvalue_n<<" "<<weighte_ilayermax<<" "<<
            //  weighte_ilayerstart<<" "<<weighte_ilayerstart_10<<" "<<weighte_ilayerend<<" "<<weighte_ilayerend_10<<endl;

        //shower weight general parameters finish	  

        // Fill Moliere radii histograms 
            mol_value = moliere(hit_energy, hit_slab, W_thicknesses, hit_x, hit_y, hit_z, hit_isMasked, masked, 0.9, shower_bool);
            radius_layer(radius90_layer_array, hit_energy, hit_slab, W_thicknesses, hit_x, hit_y, hit_z, hit_isMasked, masked, 0.9, shower_bool);
                
        // Barycenter x and y
        bary_layer(bar_layer_array, hit_energy, hit_slab, W_thicknesses, hit_x, hit_y, hit_z, hit_isMasked, masked);
        
        // Hits max distance in the same layer
        float hits_max_distance_value = hits_max_distance(hit_slab, hit_x, hit_y, hit_isMasked, masked);
        
        // Filling the tree
        b_nhit = nhit;
        b_sume = sume;
        b_weighte = weighte;
        b_bar_x = bar_xyzr[0];
        b_bar_y = bar_xyzr[1];
        b_bar_z = bar_xyzr[2];
        b_bar_r = bar_xyzr[3];
        b_mol = mol_value;
        b_MIP_Likeness = MIP_Likeness_value;
        b_hits_max_distance = hits_max_distance_value;

        b_radius90_layer_0 = radius90_layer_array[0];
            b_radius90_layer_1 = radius90_layer_array[1];
            b_radius90_layer_2 = radius90_layer_array[2];
            b_radius90_layer_3 = radius90_layer_array[3];
            b_radius90_layer_4 = radius90_layer_array[4];
            b_radius90_layer_5 = radius90_layer_array[5];
            b_radius90_layer_6 = radius90_layer_array[6];
            b_radius90_layer_7 = radius90_layer_array[7];
            b_radius90_layer_8 = radius90_layer_array[8];
            b_radius90_layer_9 = radius90_layer_array[9];
            b_radius90_layer_10 = radius90_layer_array[10];
            b_radius90_layer_11 = radius90_layer_array[11];
            b_radius90_layer_12 = radius90_layer_array[12];
            b_radius90_layer_13 = radius90_layer_array[13];
            b_radius90_layer_14 = radius90_layer_array[14];

        b_nhit_layer_0 = nhit_layer_array[0];
        b_nhit_layer_1 = nhit_layer_array[1];
        b_nhit_layer_2 = nhit_layer_array[2];
        b_nhit_layer_3 = nhit_layer_array[3];
        b_nhit_layer_4 = nhit_layer_array[4];
        b_nhit_layer_5 = nhit_layer_array[5];
        b_nhit_layer_6 = nhit_layer_array[6];
        b_nhit_layer_7 = nhit_layer_array[7];
        b_nhit_layer_8 = nhit_layer_array[8];
        b_nhit_layer_9 = nhit_layer_array[9];
        b_nhit_layer_10 = nhit_layer_array[10];
        b_nhit_layer_11 = nhit_layer_array[11];
        b_nhit_layer_12 = nhit_layer_array[12];
        b_nhit_layer_13 = nhit_layer_array[13];
        b_nhit_layer_14 = nhit_layer_array[14];

        b_nhit_layer_n_0 = nhit_layer_n_array[0];
        b_nhit_layer_n_1 = nhit_layer_n_array[1];
        b_nhit_layer_n_2 = nhit_layer_n_array[2];
        b_nhit_layer_n_3 = nhit_layer_n_array[3];
        b_nhit_layer_n_4 = nhit_layer_n_array[4];
        b_nhit_layer_n_5 = nhit_layer_n_array[5];
        b_nhit_layer_n_6 = nhit_layer_n_array[6];
        b_nhit_layer_n_7 = nhit_layer_n_array[7];
        b_nhit_layer_n_8 = nhit_layer_n_array[8];
        b_nhit_layer_n_9 = nhit_layer_n_array[9];
        b_nhit_layer_n_10 = nhit_layer_n_array[10];
        b_nhit_layer_n_11 = nhit_layer_n_array[11];
        b_nhit_layer_n_12 = nhit_layer_n_array[12];
        b_nhit_layer_n_13 = nhit_layer_n_array[13];
        b_nhit_layer_n_14 = nhit_layer_n_array[14];

        b_shower_nhit_max_layer = nhit_ilayermax;
        b_shower_nhit_start_layer = nhit_ilayerstart;
        b_shower_nhit_end_layer = nhit_ilayerend;
        b_shower_nhit_start_10_layer = nhit_ilayerstart_10;
        b_shower_nhit_end_10_layer = nhit_ilayerend_10;
        b_shower_nhit_average = nhit_shower_averagevalue;
        b_shower_nhit_max = nhit_shower_maxvalue;

        b_sume_layer_0 = sume_layer_array[0];
            b_sume_layer_1 = sume_layer_array[1];
            b_sume_layer_2 = sume_layer_array[2];
            b_sume_layer_3 = sume_layer_array[3];
            b_sume_layer_4 = sume_layer_array[4];
            b_sume_layer_5 = sume_layer_array[5];
            b_sume_layer_6 = sume_layer_array[6];
            b_sume_layer_7 = sume_layer_array[7];
            b_sume_layer_8 = sume_layer_array[8];
            b_sume_layer_9 = sume_layer_array[9];
            b_sume_layer_10 = sume_layer_array[10];
            b_sume_layer_11 = sume_layer_array[11];
            b_sume_layer_12 = sume_layer_array[12];
            b_sume_layer_13 = sume_layer_array[13];
            b_sume_layer_14 = sume_layer_array[14];

            b_sume_layer_n_0 = sume_layer_n_array[0];
            b_sume_layer_n_1 = sume_layer_n_array[1];
            b_sume_layer_n_2 = sume_layer_n_array[2];
            b_sume_layer_n_3 = sume_layer_n_array[3];
            b_sume_layer_n_4 = sume_layer_n_array[4];
            b_sume_layer_n_5 = sume_layer_n_array[5];
            b_sume_layer_n_6 = sume_layer_n_array[6];
            b_sume_layer_n_7 = sume_layer_n_array[7];
            b_sume_layer_n_8 = sume_layer_n_array[8];
            b_sume_layer_n_9 = sume_layer_n_array[9];
            b_sume_layer_n_10 = sume_layer_n_array[10];
            b_sume_layer_n_11 = sume_layer_n_array[11];
            b_sume_layer_n_12 = sume_layer_n_array[12];
            b_sume_layer_n_13 = sume_layer_n_array[13];
            b_sume_layer_n_14 = sume_layer_n_array[14];

            b_shower_sume_max_layer = sume_ilayermax;
            b_shower_sume_start_layer = sume_ilayerstart;
            b_shower_sume_end_layer = sume_ilayerend;
            b_shower_sume_start_10_layer = sume_ilayerstart_10;
            b_shower_sume_end_10_layer = sume_ilayerend_10;
            b_shower_sume_average = sume_shower_averagevalue;
            b_shower_sume_max = sume_shower_maxvalue;

        b_weighte_layer_0 = weighte_layer_array[0];
        b_weighte_layer_1 = weighte_layer_array[1];
        b_weighte_layer_2 = weighte_layer_array[2];
        b_weighte_layer_3 = weighte_layer_array[3];
        b_weighte_layer_4 = weighte_layer_array[4];
        b_weighte_layer_5 = weighte_layer_array[5];
        b_weighte_layer_6 = weighte_layer_array[6];
        b_weighte_layer_7 = weighte_layer_array[7];
        b_weighte_layer_8 = weighte_layer_array[8];
        b_weighte_layer_9 = weighte_layer_array[9];
        b_weighte_layer_10 = weighte_layer_array[10];
        b_weighte_layer_11 = weighte_layer_array[11];
        b_weighte_layer_12 = weighte_layer_array[12];
        b_weighte_layer_13 = weighte_layer_array[13];
        b_weighte_layer_14 = weighte_layer_array[14];

        b_weighte_layer_n_0 = weighte_layer_n_array[0];
        b_weighte_layer_n_1 = weighte_layer_n_array[1];
        b_weighte_layer_n_2 = weighte_layer_n_array[2];
        b_weighte_layer_n_3 = weighte_layer_n_array[3];
        b_weighte_layer_n_4 = weighte_layer_n_array[4];
        b_weighte_layer_n_5 = weighte_layer_n_array[5];
        b_weighte_layer_n_6 = weighte_layer_n_array[6];
        b_weighte_layer_n_7 = weighte_layer_n_array[7];
        b_weighte_layer_n_8 = weighte_layer_n_array[8];
        b_weighte_layer_n_9 = weighte_layer_n_array[9];
        b_weighte_layer_n_10 = weighte_layer_n_array[10];
        b_weighte_layer_n_11 = weighte_layer_n_array[11];
        b_weighte_layer_n_12 = weighte_layer_n_array[12];
        b_weighte_layer_n_13 = weighte_layer_n_array[13];
        b_weighte_layer_n_14 = weighte_layer_n_array[14];

        b_shower_weighte_max_layer = weighte_ilayermax;
        b_shower_weighte_start_layer = weighte_ilayerstart;
        b_shower_weighte_end_layer = weighte_ilayerend;
        b_shower_weighte_start_10_layer = weighte_ilayerstart_10;
        b_shower_weighte_end_10_layer = weighte_ilayerend_10;
        b_shower_weighte_average = weighte_shower_averagevalue;
        b_shower_weighte_max = weighte_shower_maxvalue;

        b_bar_x_layer_0 = bar_layer_array[0][0];
        b_bar_x_layer_1 = bar_layer_array[1][0];
        b_bar_x_layer_2 = bar_layer_array[2][0];
        b_bar_x_layer_3 = bar_layer_array[3][0];
        b_bar_x_layer_4 = bar_layer_array[4][0];
        b_bar_x_layer_5 = bar_layer_array[5][0];
        b_bar_x_layer_6 = bar_layer_array[6][0];
        b_bar_x_layer_7 = bar_layer_array[7][0];
        b_bar_x_layer_8 = bar_layer_array[8][0];
        b_bar_x_layer_9 = bar_layer_array[9][0];
        b_bar_x_layer_10 = bar_layer_array[10][0];
        b_bar_x_layer_11 = bar_layer_array[11][0];
        b_bar_x_layer_12 = bar_layer_array[12][0];
        b_bar_x_layer_13 = bar_layer_array[13][0];
        b_bar_x_layer_14 = bar_layer_array[14][0];

        b_bar_y_layer_0 = bar_layer_array[0][1];
        b_bar_y_layer_1 = bar_layer_array[1][1];
        b_bar_y_layer_2 = bar_layer_array[2][1];
        b_bar_y_layer_3 = bar_layer_array[3][1];
        b_bar_y_layer_4 = bar_layer_array[4][1];
        b_bar_y_layer_5 = bar_layer_array[5][1];
        b_bar_y_layer_6 = bar_layer_array[6][1];
        b_bar_y_layer_7 = bar_layer_array[7][1];
        b_bar_y_layer_8 = bar_layer_array[8][1];
        b_bar_y_layer_9 = bar_layer_array[9][1];
        b_bar_y_layer_10 = bar_layer_array[10][1];
        b_bar_y_layer_11 = bar_layer_array[11][1];
        b_bar_y_layer_12 = bar_layer_array[12][1];
        b_bar_y_layer_13 = bar_layer_array[13][1];
        b_bar_y_layer_14 = bar_layer_array[14][1];

        b_bar_r_layer_0 = bar_layer_array[0][2];
            b_bar_r_layer_1 = bar_layer_array[1][2];
            b_bar_r_layer_2 = bar_layer_array[2][2];
            b_bar_r_layer_3 = bar_layer_array[3][2];
            b_bar_r_layer_4 = bar_layer_array[4][2];
            b_bar_r_layer_5 = bar_layer_array[5][2];
            b_bar_r_layer_6 = bar_layer_array[6][2];
            b_bar_r_layer_7 = bar_layer_array[7][2];
            b_bar_r_layer_8 = bar_layer_array[8][2];
            b_bar_r_layer_9 = bar_layer_array[9][2];
            b_bar_r_layer_10 = bar_layer_array[10][2];
            b_bar_r_layer_11 = bar_layer_array[11][2];
            b_bar_r_layer_12 = bar_layer_array[12][2];
            b_bar_r_layer_13 = bar_layer_array[13][2];
            b_bar_r_layer_14 = bar_layer_array[14][2];

        outtree->Fill();
        
        hit_isMasked->clear();
        hit_energy->clear();
        hit_slab->clear();
        
        
    
    streamlog_out(DEBUG) << "Total energy deposit: " << totalEnergy << " GeV" <<endl;
    // Instead of filling the histograms now, we store the numbers in vectors first, then decide the binsize later
	// _evEnergyHist->Fill(totalEnergy);
	// _evHitsHist->Fill(totalHits);

}


void GetPIDParametersTreeProcessor::ShowDigitECALInfo(EVENT::LCCollection *myCollection) {
    int number = myCollection->getNumberOfElements();
    streamlog_out(DEBUG) << "TOTAL NUMBER OF HITS: " << number <<endl;
    CellIDDecoder<EVENT::CalorimeterHit> cd(myCollection);

	double totalEnergy = 0;
    int totalHits = 0;
    
    for (int i = 0; i < number; i++) {
        CalorimeterHit *ecalhit = dynamic_cast<CalorimeterHit *>(myCollection->getElementAt(i));
        
        int IJK_I = cd(ecalhit)["I"];
        int IJK_J = cd(ecalhit)["J"];
        int IJK_K = cd(ecalhit)["K"];
        float hit_energy = ecalhit->getEnergy();

        streamlog_out(DEBUG) << "\n CalorimeterHit, :" << i;
        streamlog_out(DEBUG) << " cellID-encoded=" << ecalhit->getCellID0();
        streamlog_out(DEBUG) << " I = " << IJK_I <<" mm,";
        streamlog_out(DEBUG) << " J = " << IJK_J <<" mm,";
        streamlog_out(DEBUG) << " K = " << IJK_K <<" layer,";
        streamlog_out(DEBUG) << " energy = " << hit_energy <<" GeV.\n";
		totalEnergy += hit_energy;
        totalHits++;
    }
    streamlog_out(DEBUG) << "Total energy deposit: " << totalEnergy << " GeV" <<endl;
    // Instead of filling the histograms now, we store the numbers in vectors first, then decide the binsize later
	// _evEnergyHist->Fill(totalEnergy);
	// _evHitsHist->Fill(totalHits);
	for (int i = 0; i < NUMBER_OF_LAYER; i++) {
	}
}


void GetPIDParametersTreeProcessor::processRunHeader(LCRunHeader *run)
{
}

void GetPIDParametersTreeProcessor::processEvent(LCEvent *evt) {
    try {
        streamlog_out(DEBUG) << "\n ----------------------------------------- ";
        if (_flagMcCol) {
            LCCollection *mccol = evt->getCollection(_MCColName);
            ShowMCInfo(mccol);
        }
        if (_flagEcalCol) {
            LCCollection *ecal = evt->getCollection(_ECALColName);
            ShowECALInfo(ecal);
        }
        if (_flagPixelEcalCol) {
            LCCollection *pecal = evt->getCollection(_pECALColName);
            ShowPixelECALInfo(pecal);
        }
        if (_flagDigitEcalCol) {
            LCCollection *decal = evt->getCollection(_dECALColName);
            ShowDigitECALInfo(decal);
        }
    } catch (DataNotAvailableException &e) {
        streamlog_out(DEBUG) << "Whoops!....\n";
        streamlog_out(DEBUG) << e.what();
    }
}

void GetPIDParametersTreeProcessor::check(LCEvent * evt)
{
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void GetPIDParametersTreeProcessor::end() {
    streamlog_out(MESSAGE) << "Event loop finished. Starting the fit..." <<endl;
    if (_flagMcCol) {
    }
    if (_flagEcalCol) {
    }
    
    if (_flagPixelEcalCol) {
        
    }

    if (_flagDigitEcalCol) {
    }
    streamlog_out(MESSAGE) << "Fitting done\nPlotting results..." <<endl;
}
