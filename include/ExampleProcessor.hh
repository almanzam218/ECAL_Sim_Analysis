#ifndef ExampleProcessor_h
#define ExampleProcessor_h 1
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

class ExampleProcessor : public Processor
{

public:
  virtual Processor *newProcessor() { return new ExampleProcessor; }

  ExampleProcessor();
  virtual ~ExampleProcessor();

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

//private:

  virtual void ShowMCInfo(LCCollection *col);
  virtual void ShowECALInfo(LCCollection *col);

  std::string _MCColName;
  std::string _ECALColName;


};

#endif
