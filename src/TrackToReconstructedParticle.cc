#include "TrackToReconstructedParticle.h"

TrackToReconstructedParticle aTrackToReconstructedParticle;


TrackToReconstructedParticle::TrackToReconstructedParticle() : Processor("TrackToReconstructedParticle") {

  // modify processor description
	_description = "makes reco parts from tracks given some mass/pdg info" ;


  // register steering parameters: name, description, class-variable, default value

	registerProcessorParameter( "Printing" ,
	                            "Print certain messages"  ,
	                             _printing,
	                             (int)5 ) ;

	std::vector<double> masses;
	masses.push_back(0);
	registerProcessorParameter("Masses",
				   "Mass of the PDG code assigned",
				    _masses,
				    masses);

	std::vector<int> charges;
	charges.push_back(0);
	registerProcessorParameter("Charges",
				   "Charge of the PDG code assigned",
				    _charges,
				    charges);

	std::vector<int> pdgs;
	pdgs.push_back(0);
	registerProcessorParameter("Pdgs",
				   "PDG codes to be assigned to tracks",
				    _pdgs,
				    pdgs);

   	std::string inputTrackCollectionName = "MarlinTrkTracks";
  	registerInputCollection( LCIO::TRACK,
                                 "InputTrackCollectionName" ,
                                 "Input Track Collection Name " ,
                                 _inputTrackCollectionName,
                                 inputTrackCollectionName);  

  	std::string outputParticleCollectionName = "NewPfoCol";
  	registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                             	"OutputParticleCollectionName" ,
			     	"Output Particle Collection Name "  ,
                             	_outputParticleCollectionName,
                             	outputParticleCollectionName);

}

void TrackToReconstructedParticle::init() {

  streamlog_out(DEBUG) << "   init called  "
                       << std::endl ;


  // usually a good idea to
  printParameters() ;
  nEvt = 0;
 
  BField = marlin::Global::GEAR->getBField().at(gear::Vector3D(0.,0.,0.)).z();
//  gROOT->ProcessLine("#include <vector>");
}

void TrackToReconstructedParticle::processRunHeader( LCRunHeader* run) {
  streamlog_out(MESSAGE) << " processRunHeader "  << run->getRunNumber() << std::endl ;
}
void TrackToReconstructedParticle::printTrack(Track* t){
	std::cout<<"Track: (d0,phi,ome,z0,tanL) "<< 
		t->getD0()<<" "<<
		t->getPhi()<<" "<<
		t->getOmega()<<" "<<
		t->getZ0()<<" "<<
		t->getTanLambda()<<std::endl;
}
void TrackToReconstructedParticle::printReconstructedParticle(ReconstructedParticle* p){
	const double* mom = p->getMomentum();	
	std::cout<<"Particle "<< p->getType() <<": "<<
	"(Px,Py,Pz,E,M,q) "<<
	mom[0]<< " "<<mom[1]<< " "<<mom[2]<< " "
	<<p->getEnergy()<<" "<<p->getMass()<<" "
	<<p->getCharge()<<std::endl;
		
	
}
bool TrackToReconstructedParticle::FindTracks( LCEvent* evt ) {

  bool tf = false;

  // clear old vector
  _trackvec.clear();
  typedef const std::vector<std::string> StringVec ;
  StringVec* strVec = evt->getCollectionNames() ;
  for(StringVec::const_iterator itname=strVec->begin(); itname!=strVec->end(); itname++){
    if(*itname==_inputTrackCollectionName){
      LCCollection* col = evt->getCollection(*itname);
      unsigned int nelem = col->getNumberOfElements();
      tf = true;
      for(unsigned int i=0;i<nelem;i++){
        Track* track = dynamic_cast<Track*>(col->getElementAt(i));
        _trackvec.push_back(track);
      }
    }
  }

  if(_printing>1)std::cout << "FindTracks : " << tf << std::endl;

  return tf;
}
std::vector<double> TrackToReconstructedParticle::getTrackPxPyPz(Track* t){
	const double c = 2.99792458e8; // m*s^-1        
  	const double mm2m = 1e-3;
  	const double eV2GeV = 1e-9;
  	const double eB = BField*c*mm2m*eV2GeV;
 
 	
	double cosLambda = 1 / sqrt(1 + t->getTanLambda()*t->getTanLambda() );
	double P = (eB/fabs(t->getOmega()))/cosLambda;
	double sinLambda = t->getTanLambda()*cosLambda;
	double cosPhi = cos(t->getPhi());
	double sinPhi = sin(t->getPhi());
	double px = P*cosLambda*cosPhi;
	double py = P*cosLambda*sinPhi;
	double pz = P*sinLambda;
	std::vector<double> txtytz;
	txtytz.push_back(px);
	txtytz.push_back(py);
	txtytz.push_back(pz);
	return txtytz;
}

ReconstructedParticle* TrackToReconstructedParticle::makePartFromTrack(Track* t, double mass, int charge, int pdg){

		ReconstructedParticleImpl* p = new ReconstructedParticleImpl();
		ParticleIDImpl* newPDG = new ParticleIDImpl();
		newPDG->setPDG(pdg);
		newPDG->setLikelihood(1.0);
		
		//get pxpypz from track, calculate E from TLorentzVector
		std::vector<double> pvec = getTrackPxPyPz(t);
	
		TLorentzVector v;
		v.SetXYZM(pvec.at(0),pvec.at(1),pvec.at(2),mass );
		
		float* mom = new float[3];

		mom[0] = v.Px();
		mom[1] = v.Py();
 		mom[2] = v.Pz();	
		
		p->setMomentum(mom);
		p->setEnergy( v.E() );

		p->setMass(mass);
		p->setCharge(charge);
		p->addParticleID(newPDG);
		p->setParticleIDUsed(newPDG);
		p->setType(pdg);
		//dont worry about setting the cov in the 
		//reconstructedparticle, just only use the
		//track covariance matrix
		return p;	
}

void TrackToReconstructedParticle::processEvent( LCEvent * evt ) {
 //FindMCParticles(evt);
// = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
 FindTracks(evt);
  
  LCCollectionVec * partCollection = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);

  streamlog_out(MESSAGE) << " start processing event " << std::endl;
 

 	std::cout<<"track vector size "<<_trackvec.size()<<std::endl;

	//loop over tracks
	ReconstructedParticle* part;
	for(unsigned int i =0; i<_trackvec.size(); i++){
		//for each mass/pdg assignment create a particle
		//make sure the charge indicated is consistent with the track charge
		for( unsigned int j=0; j< _masses.size(); i++){
			//if the charges match continue
			if( _trackvec.at(i)->getOmega() * _charges.at(j) > 0){
				part = makePartFromTrack(_trackvec.at(i), masses.at(j) , charges.at(j), pdgs.at(j) );					
				partCollection->addElement( part );
				std::cout<<"Particle Created:"<<std::endl;
				printTrack(_trackvec.at(i));
				printReconstructedParticle( part );
			}
		}
				

	}

  

  nEvt++;

  // Add new collection to event

  evt->addCollection(partCollection , _outputParticleCollectionName.c_str() ); 
 std::cout << "======================================== event " << nEvt << std::endl ;

}
void TrackToReconstructedParticle::end(){

}

