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

// ----- include for verbosity dependent logging ---------
#include "marlin/VerbosityLevels.h"
#include "marlin/StringParameters.h"
#define SLM streamlog_out(MESSAGE)

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

class GetMIPProcessor : public Processor
{
private:
    const static int NUMBER_OF_LAYER = 15;
    const static int NUMBER_OF_CELLX = 64;
    const static int NUMBER_OF_CELLY = 32;
    float FIT_INIT_AMP = 5000;//events
    float FIT_INIT_MPV = 0.1e-3;//GeV
    float FIT_INIT_SIG = 0.01e-3;//GeV
    float FIT_RANGE_MIN = 0;
    float FIT_RANGE_MAX = 1e-3;//GeV

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
    double _layerFitParams[NUMBER_OF_LAYER][3];
    double _pixelLayerFitParams[NUMBER_OF_LAYER][3];
  
    TF1* landauFunc = new TF1("fitLandauFunc", "landau", FIT_RANGE_MIN, FIT_RANGE_MAX);
    TH1* _fittedMIP;
    TH1* _fittedPixelMIP;

private:
  virtual void ShowMCInfo(LCCollection *col);
  virtual void ShowECALInfo(LCCollection *col);
    virtual void ShowPixelECALInfo(LCCollection *col);

  std::string _MCColName;
  std::string _ECALColName;
    std::string _pECALColName;

};

#endif
