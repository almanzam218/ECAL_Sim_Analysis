#ifndef GetPIDParametersTreeProcessor_h
#define GetPIDParametersTreeProcessor_h 1
#include <iomanip>
#include <EVENT/LCRelation.h>
#include "marlin/Processor.h"
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/TrackStateImpl.h>
#include <EVENT/Track.h>
#include <EVENT/Vertex.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <UTIL/LCRelationNavigator.h>
#include <UTIL/PIDHandler.h>
#include <UTIL/CellIDDecoder.h>
#include <cmath>
#include <string>
#include <vector>
#include <array>

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TImage.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include <TF1.h>
#include <TRandom.h>
#include <TFile.h>
#include <TTree.h>
// #include <TLorentzVector.h>
//MAS: Included this three lines
#ifdef MARLIN_AIDA //AIDA
#include <marlin/AIDAProcessor.h>
#endif

using namespace lcio;
using namespace marlin;
using namespace std;

class GetPIDParametersTreeProcessor : public Processor {

private:
    const static int NUMBER_OF_LAYER = 15;
    const static int NUMBER_OF_CELLX = 64;
    const static int NUMBER_OF_CELLY = 32;
    float FIT_INIT_AMP = 5000;//events
    float FIT_INIT_MPV = 0.1e-3;//GeV
    float FIT_INIT_SIG = 0.01e-3;//GeV
    float E_RANGE_MIN = 0;
    float E_RANGE_MAX = 0.05;//GeV

public:
  virtual Processor *newProcessor() { return new GetPIDParametersTreeProcessor; }

  GetPIDParametersTreeProcessor();
  virtual ~GetPIDParametersTreeProcessor();

  /** Called at the begin of the job before anything is read.
   */
  virtual void init();

  /** Called for every run.
   */
  virtual void processRunHeader(LCRunHeader *run);

  /** Called for every event - the working horse.
   */
  virtual void processEvent(LCEvent *evt);

  virtual void check(LCEvent *evt);

  virtual void end();

  // Histogram definitions for ECALHit class

    // MC particle
    // Monolithic calorimeter hits
    // Pixelised calorimeter hits
    // Digitised calorimeter hits
    // double _layerFitParams[NUMBER_OF_LAYER][4];
    // double energyRes;
    

private:
  virtual void ShowMCInfo(LCCollection *col);
  virtual void ShowECALInfo(LCCollection *col);
    virtual void ShowPixelECALInfo(LCCollection *col);
    virtual void ShowDigitECALInfo(LCCollection *col);

  std::string _MCColName;
  std::string _ECALColName;
    std::string _pECALColName;
    std::string _dECALColName;
    
    bool _flagMcCol = false;
    bool _flagEcalCol = false;
    bool _flagPixelEcalCol = false;
    bool _flagDigitEcalCol = false;
    float runEnergy = -1;
    Float_t b_ecal_interaction;
    Float_t b_nhit, b_sume, b_weighte, b_bar_x, b_bar_y, b_bar_z, b_bar_r;
    Float_t b_mol;
    Float_t b_MIP_Likeness;
    Float_t b_hits_max_distance;
    Float_t b_radius90_layer_0, b_radius90_layer_1, b_radius90_layer_2, b_radius90_layer_3, b_radius90_layer_4, b_radius90_layer_5, b_radius90_layer_6, b_radius90_layer_7, b_radius90_layer_8, b_radius90_layer_9, b_radius90_layer_10, b_radius90_layer_11, b_radius90_layer_12, b_radius90_layer_13, b_radius90_layer_14;
    Float_t b_bar_x_layer_0, b_bar_x_layer_1, b_bar_x_layer_2, b_bar_x_layer_3, b_bar_x_layer_4, b_bar_x_layer_5, b_bar_x_layer_6, b_bar_x_layer_7, b_bar_x_layer_8, b_bar_x_layer_9, b_bar_x_layer_10, b_bar_x_layer_11, b_bar_x_layer_12, b_bar_x_layer_13, b_bar_x_layer_14;
    Float_t b_bar_y_layer_0, b_bar_y_layer_1, b_bar_y_layer_2, b_bar_y_layer_3, b_bar_y_layer_4, b_bar_y_layer_5, b_bar_y_layer_6, b_bar_y_layer_7, b_bar_y_layer_8, b_bar_y_layer_9, b_bar_y_layer_10, b_bar_y_layer_11, b_bar_y_layer_12, b_bar_y_layer_13, b_bar_y_layer_14;
    Float_t b_bar_r_layer_0, b_bar_r_layer_1, b_bar_r_layer_2, b_bar_r_layer_3, b_bar_r_layer_4, b_bar_r_layer_5, b_bar_r_layer_6, b_bar_r_layer_7, b_bar_r_layer_8, b_bar_r_layer_9, b_bar_r_layer_10, b_bar_r_layer_11, b_bar_r_layer_12, b_bar_r_layer_13, b_bar_r_layer_14;
    Float_t b_shower_nhit_max_layer, b_shower_nhit_start_layer, b_shower_nhit_end_layer, b_shower_nhit_start_10_layer, b_shower_nhit_end_10_layer, b_shower_nhit_average, b_shower_nhit_max;
    Float_t b_shower_sume_max_layer, b_shower_sume_start_layer, b_shower_sume_end_layer, b_shower_sume_start_10_layer, b_shower_sume_end_10_layer, b_shower_sume_average, b_shower_sume_max;
    Float_t b_shower_weighte_max_layer, b_shower_weighte_start_layer, b_shower_weighte_end_layer, b_shower_weighte_start_10_layer, b_shower_weighte_end_10_layer, b_shower_weighte_average, b_shower_weighte_max;
    Float_t b_nhit_layer_0, b_nhit_layer_1, b_nhit_layer_2, b_nhit_layer_3, b_nhit_layer_4, b_nhit_layer_5, b_nhit_layer_6, b_nhit_layer_7, b_nhit_layer_8, b_nhit_layer_9, b_nhit_layer_10, b_nhit_layer_11, b_nhit_layer_12, b_nhit_layer_13, b_nhit_layer_14;
    Float_t b_nhit_layer_n_0, b_nhit_layer_n_1, b_nhit_layer_n_2, b_nhit_layer_n_3, b_nhit_layer_n_4, b_nhit_layer_n_5, b_nhit_layer_n_6, b_nhit_layer_n_7, b_nhit_layer_n_8, b_nhit_layer_n_9, b_nhit_layer_n_10, b_nhit_layer_n_11, b_nhit_layer_n_12, b_nhit_layer_n_13, b_nhit_layer_n_14;
    Float_t b_weighte_layer_0, b_weighte_layer_1, b_weighte_layer_2, b_weighte_layer_3, b_weighte_layer_4, b_weighte_layer_5, b_weighte_layer_6, b_weighte_layer_7, b_weighte_layer_8, b_weighte_layer_9, b_weighte_layer_10, b_weighte_layer_11, b_weighte_layer_12, b_weighte_layer_13, b_weighte_layer_14;
    Float_t b_weighte_layer_n_0, b_weighte_layer_n_1, b_weighte_layer_n_2, b_weighte_layer_n_3, b_weighte_layer_n_4, b_weighte_layer_n_5, b_weighte_layer_n_6, b_weighte_layer_n_7, b_weighte_layer_n_8, b_weighte_layer_n_9, b_weighte_layer_n_10, b_weighte_layer_n_11, b_weighte_layer_n_12, b_weighte_layer_n_13, b_weighte_layer_n_14;    
    Float_t b_sume_layer_0, b_sume_layer_1, b_sume_layer_2, b_sume_layer_3, b_sume_layer_4, b_sume_layer_5, b_sume_layer_6, b_sume_layer_7, b_sume_layer_8, b_sume_layer_9, b_sume_layer_10, b_sume_layer_11, b_sume_layer_12, b_sume_layer_13, b_sume_layer_14;
    Float_t b_sume_layer_n_0, b_sume_layer_n_1, b_sume_layer_n_2, b_sume_layer_n_3, b_sume_layer_n_4, b_sume_layer_n_5, b_sume_layer_n_6, b_sume_layer_n_7, b_sume_layer_n_8, b_sume_layer_n_9, b_sume_layer_n_10, b_sume_layer_n_11, b_sume_layer_n_12, b_sume_layer_n_13, b_sume_layer_n_14;


};

#endif
