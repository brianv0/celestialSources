#include "astro/SkyDir.h"
#include "astro/GPS.h"
#include "FluxSvc/IFluxSvc.h"


#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"

#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/MonteCarlo/McParticle.h"
#include "Event/MonteCarlo/McIntegratingHit.h"
#include "Event/MonteCarlo/McPositionHit.h"

#include "Event/Recon/TkrRecon/TkrFitTrack.h" //aggiunto 
#include "Event/Recon/TkrRecon/TkrKalFitTrack.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"
#include "Event/Recon/TkrRecon/TkrClusterCol.h" //max
#include "Event/Recon/TkrRecon/TkrCluster.h" //max

#include "Event/Recon/CalRecon/CalCluster.h"
//include files for ROOT...
#include "TTree.h"
#include "TH1D.h"
#include "TObjArray.h"
#include "TFile.h"
#include "TMath.h"

//#include "TkrRecon/Cluster/TkrMakeClusters.h" //max
#include "GaudiKernel/DataSvc.h" //max
//#include <map>
#include <vector>
#include <fstream>
#include <sstream>
#include "CLHEP/Vector/Rotation.h"

#include "LatGRBAlert/PhotonInfo.h"
#include "LatGRBAlert/ClusterData.h"
#include "LatGRBAlert/LatGRBAlertConstants.h"

bool SaveFullRecon   = true;
int GRBAlertFlag     = 1; //0 = none, 1 = FR, 2 = MC

struct myEvent
{  
  double time;
  double energy;
  double theta;
  double phi;
  double l;
  double b;
};

struct GRBAlertData
{ 
  double tstart;
  double tend;
  double like;
};

const double deg   =  180.0 / M_PI;

astro::SkyDir GalacticDirection(Hep3Vector GlastDirection, double time)
{
  HepRotation GlastToGalactic = GPS::instance()->transformCelToGlast(time).inverse();
  astro::SkyDir GalDir(GlastToGalactic * GlastDirection);
  return GalDir;
}

class TDSReadFluxAlg : public Algorithm
{	
public:
  //! A constructor of this type has to be provide
  TDSReadFluxAlg(const std::string& name, ISvcLocator* pSvcLocator);
  
  /*! \brief Initializes the Algorithm
    
    It sets the correct pointers and builds a root tree for saving the data
  */
  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();   
  
private:
  HepRotation Glast2Gal;  
  Hep3Vector dir,GalDirVec;
    
  StatusCode readMCEvent();
  StatusCode readEvent();
  StatusCode GRBAlert(myEvent);
  myEvent MCEvent,FREvent;  

  int NumReconTracks;
  

  std::vector<int> Hits;
  TTree *FRevents,*MCevents;

  IFluxSvc* m_fsvc; // pointer to the flux Service
  
  //  TH1D *hitmap[numPlanes];//,*hitmap0y,*hitmap1x,*hitmap1y,*hitmap2x,*hitmap2y;
  int numPlanes,numTowers;
  bool triggerFlag;
  std::vector<std::string> m_save_file;
  //  std::string m_nameoutfile;
  //std::string m_save_opt; 
  
  PhotonInfo a_photon;
  std::vector<PhotonInfo> PhotonsData;
  std::vector<GRBAlertData> JointLike;
  int count_buff;
  //////////////////////////////////////////////////
  ClusterData *clusterData;
  //////////////////////////////////////////////////
};


static const AlgFactory<TDSReadFluxAlg>  Factory;
const IAlgFactory& TDSReadFluxAlgFactory = Factory;


TDSReadFluxAlg::TDSReadFluxAlg(const std::string& name, 
			       ISvcLocator* pSvcLocator) : 
  Algorithm(name, pSvcLocator)
{
  declareProperty("savefile", m_save_file);
  //  declareProperty("nameoutput", m_nameoutfile);
 }

StatusCode TDSReadFluxAlg::initialize()
{
  NumReconTracks=0;
  StatusCode sc = StatusCode::SUCCESS;
  MsgStream log(msgSvc(), name());
  setProperties();
  sc = service("FluxSvc", m_fsvc);
  //////////////////////////////////////////////////
  /*
    IService*   iService = 0; // questa parte serve per i cluster...
    sc        = serviceLocator()->getService("TkrGeometrySvc", iService, true);  
    sc        = serviceLocator()->getService("EventDataSvc", iService);
    m_dataSvc = dynamic_cast<DataSvc*>(iService);  
    //Geometria del TKR..
    sc = service("TkrGeometrySvc", m_TkrGeo, true);
    if (sc.isFailure()) 
    {
    log << MSG::ERROR << "TkrGeometrySvc is required for this algorithm." 
    << endreq;
    return sc;
    }
    
    m_Alignment = m_TkrGeo->getTkrAlignmentSvc();
    numTowers = m_TkrGeo->numXTowers()*m_TkrGeo->numYTowers();
    numPlanes =  2*m_TkrGeo->numPlanes()*numTowers;
  */
  
  MCevents = new TTree("MCEvents","MCEvents");    
  MCevents->Branch("time",&MCEvent.time,"time/D");  
  MCevents->Branch("energy",&MCEvent.energy,"energy/D");  
  MCevents->Branch("theta",&MCEvent.theta,"theta/D");  
  MCevents->Branch("phi",&MCEvent.phi,"phi/D");
  MCevents->Branch("l",&MCEvent.l,"l/D");
  MCevents->Branch("b",&MCEvent.b,"b/D");
  
  FRevents = new TTree("FREvents","FREvents");    
  FRevents->Branch("time",&FREvent.time,"time/D");  
  FRevents->Branch("energy",&FREvent.energy,"energy/D");  
  FRevents->Branch("theta",&FREvent.theta,"theta/D");  
  FRevents->Branch("phi",&FREvent.phi,"phi/D");
  FRevents->Branch("l",&FREvent.l,"l/D");
  FRevents->Branch("b",&FREvent.b,"b/D");
  
  /*
    MCevents->Branch("MCevents",&MCEvent,
    "time/D:energy:Theta:Phi:l:b");
    
    FRevents = new TTree("FREvents","FREvents");  
    FRevents->Branch("FRevents",&FREvent,
    "time/D:energy:Theta:Phi:l:b");
  */
  if(GRBAlertFlag>0) clusterData = new ClusterData;
  //////////////////////////////////////////////////
  return sc;
}

StatusCode TDSReadFluxAlg::execute()
{
  MsgStream log(msgSvc(), name());
  StatusCode sc = StatusCode::SUCCESS;
  sc = readMCEvent();
  sc = readEvent();
  return sc;
}

StatusCode TDSReadFluxAlg::readMCEvent()
{
  StatusCode sc = StatusCode::SUCCESS;

  MsgStream log(msgSvc(), name());
  // ---> MC PARTICLES
  SmartDataPtr<Event::EventHeader> evt(eventSvc(), EventModel::EventHeader);
  if(!evt) return sc;
  //Get the time
  double time = (double) evt->time();
  MCEvent.time= time;
  
  SmartDataPtr<Event::McParticleCol> 
    particles(eventSvc(), EventModel::MC::McParticleCol);
  

  if (particles) 
    {
      Event::McParticleCol::const_iterator p;
      p = particles->begin(); // This is the mother...
      
      MCEvent.energy = (*p)->finalFourMomentum().e();
      dir= -1*(*p)->finalFourMomentum().vect().unit();
      
      MCEvent.theta     =  acos(dir.cosTheta())*deg; // note the (-)
      MCEvent.phi       =  dir.phi()*deg; //conversione in phi-Xml
      
      MCEvent.l = GalacticDirection(dir,time).l();
      MCEvent.b = GalacticDirection(dir,time).b();
      
      //////////////////////////////////////////////////
      log << MSG::INFO << "Retrieved McParticles from the TDS at time = " << MCEvent.time <<" s "<< endreq;
      log << MSG::INFO <<" MC energy           = "<<MCEvent.energy<<" (MeV)"<< endreq;
      log << MSG::INFO <<" MC Direction: Theta = "<<MCEvent.theta<<" Phi = "<<MCEvent.phi<< endreq;
      log << MSG::INFO <<" MC Direction: l     = "<<MCEvent.l<<" b = "<<MCEvent.b<< endreq;
      if(GRBAlertFlag == 2) sc = GRBAlert(MCEvent);      
      if(MCEvent.energy>0) MCevents->Fill();
    }  
  return  StatusCode::SUCCESS;
}

  
StatusCode TDSReadFluxAlg::readEvent()
{
  MsgStream log(msgSvc(), name());
  StatusCode sc = StatusCode::SUCCESS;
  // ---> FULL RECONSTRUCTION
  // ---> 1) Calorimeter:
  SmartDataPtr<Event::CalClusterCol> 
    calClusters(eventSvc(),  EventModel::CalRecon::CalClusterCol);
  
  if(calClusters) // look at the recon energy
    {
      FREvent.time = MCEvent.time; 
      FREvent.energy = calClusters->getCluster(0)->getEnergyCorrected(); 
      log << MSG::INFO <<" FR Energy           = "<<FREvent.energy<<" (MeV)"<< endreq;
    }
  
  // ---> 2) Tracker:
  
  SmartDataPtr<Event::TkrClusterCol> 
    tkrClusters(eventSvc(),  EventModel::TkrRecon::TkrClusterCol);
  
  SmartDataPtr<Event::TkrVertexCol>  
    tracks(eventSvc(),EventModel::TkrRecon::TkrVertexCol);
  
  
  if(tkrClusters)	
    {
      if (tracks->size() != 0)
	{
	  const Event::TkrVertex *track = (tracks->front());
	  
	  dir = -1*track->getDirection();
	  
	  FREvent.theta = acos(dir.cosTheta())*deg; // note the (-)
	  FREvent.phi   = dir.phi()*deg; //conversione in Phi-Xml
	  FREvent.l = GalacticDirection(dir,FREvent.time).l();
	  FREvent.b = GalacticDirection(dir,FREvent.time).b();
	  
	  log << MSG::INFO <<" FR Direction: Theta = "<<FREvent.theta<<" Phi = "<<FREvent.phi<< endreq;
	  log << MSG::INFO <<" FR Direction: l     = "<<FREvent.l<<" b = "<<FREvent.b<< endreq;
	  
	  if(GRBAlertFlag == 1)  sc = GRBAlert(FREvent);
	  NumReconTracks++;
	  std::cout<<" Number of reconstructed tracks = "<<NumReconTracks<<std::endl;
	  FRevents->Fill(); 
	}
    }
  
  if(triggerFlag) FRevents->Fill();       
  // FINE FULL RECON
  //  OBrecon *myOB = new OBrecon();
  
  return sc;
}


StatusCode TDSReadFluxAlg::finalize()
{
  std::cout<<"--------------------------------------------------"<<std::endl;
  std::cout<<" Number of reconstructed tracks = "<<NumReconTracks<<std::endl;
  std::cout<<"--------------------------------------------------"<<std::endl;
  
  if(JointLike.size()>0)
    {
      std::ofstream os("trigger_Likelihood.lis");
      for (std::vector<GRBAlertData>::iterator it = JointLike.begin(); 
	   it!=JointLike.end(); ++it)
	{
	  os << it->tstart<<"  "<<it->tend<<"  "<<it->like<<std::endl; 
	  std::cout<<it->tstart<<"  "<<it->tend<<"  "<<it->like<<std::endl; 
	}
    }
  
  if(m_save_file.empty()) return StatusCode::SUCCESS;
  std::vector<std::string>::iterator itr;
  for(itr=m_save_file.begin();itr != m_save_file.end();++itr)
    {
      if((*itr)=="root") 
	{
	  std::cout<<"Saving Events.root ..."<<std::endl;
	  TFile *f = new TFile("Events.root","RECREATE"); //real events
	  MCevents->Write();  
	  FRevents->Write();
	  f->Close();
	  std::cout<<"Events.root Saved !"<<std::endl;
	}
      // Write the ascii file
      /*
	if((*itr)=="ascii") 
	{
	long sz = m_data.size();
	if (sz > 0)   
	{
	std::ofstream os("grb_detected.lis");
	os << "nEvents: " << sz << std::endl;
	#ifdef WIN32
	for (std::vector<Data>::iterator it=m_data.begin(); 
	it!=m_data.end(); ++it)
	os << std::setw(12) <<
	std::fixed << 
	it->time << "  " << 
	it->energy << "  " <<
	it->theta << "  " << 
	it->phi << std::endl;
	#else
	for (std::vector<Data>::iterator it=m_data.begin(); 
	it!=m_data.end(); ++it)
	os << std::setw(12) << 
	it->time << "  " << 
	it->energy << "  " << 
	it->theta << "  " << 
	it->phi << std::endl;
	#endif
	}
	}
      */
    }
  return StatusCode::SUCCESS;
}

StatusCode TDSReadFluxAlg::GRBAlert(myEvent Evt)
{
  //////////////////////////////////////////////////
  ///           LAT GRB ALERT
  //////////////////////////////////////////////////
  PhotonInfo myPhoton;
 
  std::cout<<Evt.time<<" "<<Evt.energy<<std::endl;

  myPhoton.setTime(Evt.time);
  myPhoton.setEnergy(Evt.energy);
  myPhoton.setTheta(Evt.l);
  myPhoton.setPhi(Evt.b);
  
  PhotonsData.push_back(myPhoton);    
  GRBAlertData AlertData;
  
  std::cout<<PhotonsData[0].time()<<" "<<PhotonsData.back().time()<<std::endl;

  if((int) PhotonsData.size() == lattrigcst::nrange)
    {
      double tstart = PhotonsData[0].time();
      double tend   = PhotonsData.back().time();
      double like   = clusterData->triggerLikelihood(PhotonsData);
      std::cout<<tstart<<" "<<tend<<" "<<like<<std::endl;
      if(like>100)
	{
	  std::cout << " *********************************************** " << std::endl;
	  std::cout << " ********** GRB TRIGGER ON  ("<<like
		    <<" ) TIME = "<<tstart<<" ********** " << std::endl;
	  std::cout << "************************************************ " << std::endl;
	}
      
      AlertData.tstart = tstart;
      AlertData.tend   = tend;
      AlertData.like   = like;

      JointLike.push_back(AlertData);

      std::cout << " --> Likelihood  ("<<tstart<<","<<tend<<") = " <<like << std::endl;
      std::copy(&PhotonsData[lattrigcst::nmove], 
		&PhotonsData[lattrigcst::nmove+lattrigcst::nrange], PhotonsData.begin());
      PhotonsData.resize(lattrigcst::nrange-lattrigcst::nmove);
    }
  std::cout << " Trigger buffer = " << PhotonsData.size()<<" / "<<lattrigcst::nrange
	    << std::endl;
  
  return StatusCode::SUCCESS;
}
