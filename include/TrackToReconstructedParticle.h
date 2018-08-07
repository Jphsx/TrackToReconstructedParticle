#include "marlin/Processor.h"
#include "EVENT/Track.h"
#include "lcio.h"
#include "TFile.h"
#include <vector>
#include "IMPL/LCCollectionVec.h"
#include "IMPL/ParticleIDImpl.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include <marlin/Global.h>
#include "gear/BField.h"
#include "TLorentzVector.h"

using namespace lcio;

	/** TrackToReconstructedParticle:<br>
 *
 * 
 * @author Justin Anguiano, University of Kansas
 * 
 */

 class TrackToReconstructedParticle : public marlin::Processor {

 public:

 virtual marlin::Processor*  newProcessor() { return new TrackToReconstructedParticle ; }

  TrackToReconstructedParticle(const TrackToReconstructedParticle&) = delete ;
  TrackToReconstructedParticle& operator=(const TrackToReconstructedParticle&) = delete ;

  TrackToReconstructedParticle() ;

  /** Called at the beginning of the job before anything is read.
   *  Use to initialize the proscessor, e.g. book histograms.
   */
  virtual void init() ;
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;

  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ;


  /** Called after data processing for clean up.
   */
  virtual void end() ;

  bool FindTracks(LCEvent* evt);
 
  std::vector<double> getTrackPxPyPz(Track* t);
  ReconstructedParticle* makePartFromTrack(Track* t, double mass, int charge, int pdg);
  
  //printing utility
  void printTrack(Track* t);
  void printReconstructedParticle(ReconstructedParticle* p);

  protected:
  int nEvt{};
  
  //vector to hold the tracks for the event
  std::vector<Track*> _trackvec{};
  int   _printing{};

   //need to no BField to calculate stuff
   double BField{};

  // particle info to be assigned to tracks
   std::vector<double> _masses{};
   std::vector<int> _charges{};
   std::vector<int> _pdgs{};
 
  
// _inputTrackCollectionName 
  std::string _outputParticleCollectionName{};
  std::string _inputTrackCollectionName{};
//  std::string m_rootFile{};
};
