#ifndef GetMIPProcessor_h
#define GetMIPProcessor_h 1
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

class GetMIPProcessor : public Processor
{
private:
    int NUMBER_OF_LAYER = 15;
    int NUMBER_OF_CELLX = 64;
    int NUMBER_OF_CELLY = 32;

public:
  virtual Processor *newProcessor() { return new GetMIPProcessor; }

  GetMIPProcessor();
  virtual ~GetMIPProcessor();

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

  TH1* _xHist;
  TH1* _yHist;
  TH1* _zHist;
  TH2* _xyHist;
  TH1* _cellEnergyHist;
  TH1F* _energyInLayerSi[NUMBER_OF_LAYER];
  TH1F* _energyInPixelLayerSi[NUMBER_OF_LAYER];
  double _layerFitParams[NUMBER_OF_LAYER][4];
  double _pixelLayerFitParams[NUMBER_OF_LAYER][4];

private:
  virtual void ShowMCInfo(LCCollection *col);
  virtual void ShowECALInfo(LCCollection *col);
    virtual void ShowPixelECALInfo(LCCollection *col);

  std::string _MCColName;
  std::string _ECALColName;
    std::string _pECALColName;

};

#endif
