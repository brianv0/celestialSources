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
#include "Event/Recon/TkrRecon/TkrVertex.h"
#include "Event/Recon/CalRecon/CalCluster.h"

#include "Event/Digi/CalDigi.h"
#include "idents/CalXtalId.h"
//include files for ROOT...
#include "TTree.h"
#include "TBranch.h"
#include "TObjArray.h"
#include "TFile.h"
// Include files
#include "FluxSvc/IFluxSvc.h"
//#include <map>
#include <vector>
#include <fstream>
#include "CLHEP/Vector/Rotation.h"
/*! @class TDSReadFluxAlg
  @brief Takes the data relatives to the incoming particle from the TDS
  and write the output in files.
 
  This method accesses to the transient data store to retrive information
  about the incoming particles (mother particles) and about their trigger 
  status and their reconstruction.  
  From the TDS:
  - The time (from EventHeader)
  - The "real" energy (from McParticleCol)
  - The reconstructed energy (from CalClusterCol)
  - The "real" direction (from McParticleCol)
  - The reconstructed direction (from TkrVertexCol)
  Then it compute the galactic direction (using \e transformGlastToGalactic(time)
  from \e fluxSvc.
  Finally this algorithm save the "real" data and the "recon" data in 
  output files.
  From the jobOptions file it is possible to select the format: ROOT tree, 
  ascii file, or both. 
  To test the algorithm, run: 
  \verbatim
  test_GRB.exe ../src/test/TDSreadFluxOptions.txt
  \endverbatim
  
  @author Nicola Omodei, Francesco Longo, Sandhia Bansal
*/

const double deg = 180.0 / M_PI;

/*! @struct Data 
  \brief This structure stores the data
  
  The data that we read from the TDS and we whant to store are:
  \param time the time at which the particle arrives to the satellite
  \param energy is the energy of the particle
  \param theta relative to Glast
  \param phi relative to glast
  
  @author Nicola Omodei, Francesco Longo, Sandhia Bansal
*/
struct Data
{
	float  time;
	float  energy;
	float  theta;
	float  phi;
};

std::pair<double,double> GalDir(HepVector3D GalDirVec)
{
	double x_g=GalDirVec.x();
	double y_g=GalDirVec.y();
	double z_g=GalDirVec.z();

	if(pow(x_g,2)<=1.e-15) x_g=0.0;
	if(pow(y_g,2)<=1.e-15) y_g=0.0;
	if(pow(z_g,2)<=1.e-15) z_g=0.0;

	double b = 360./M_2PI*asin(y_g);
	double l = 360./M_2PI*atan2(x_g,z_g);
	if(pow(b,2)<1.0e-15) b=0.;
	if(pow(l,2)<1.0e-15) l=0.;

	return std::make_pair(l,b);

}

class TDSReadFluxAlg : public Algorithm
{	
public:
    TDSReadFluxAlg(const std::string& name, ISvcLocator* pSvcLocator);
    
    StatusCode initialize();
   
    StatusCode execute();
    
    StatusCode finalize();   
        
private:
    IFluxSvc* m_fsvc; // pointer to the flux Service
    HepRotation Glast2Gal;  
    HepVector3D dir,GalDirVec;
    StatusCode readTDSData();

    double energy,time,cos_theta,phi,l,b;
    double energy_recon,time_recon,cos_theta_recon,phi_recon,l_recon,b_recon,rate_recon;
    int trigger_flag;
    std::vector<Data> m_data;
    
    TTree *events;
    TTree *events_recon;
   
    TFile *f;
    TFile *f_recon;
    int counts;
    std::vector<std::string> m_save_file;
//  std::string m_save_opt;
};

  
static const AlgFactory<TDSReadFluxAlg>  Factory;
const IAlgFactory& TDSReadFluxAlgFactory = Factory;


TDSReadFluxAlg::TDSReadFluxAlg(const std::string& name, ISvcLocator* pSvcLocator) : 
Algorithm(name, pSvcLocator)
{
	declareProperty("savefile", m_save_file);
//	declareProperty("saveoption", m_save_opt);

}

StatusCode TDSReadFluxAlg::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    
    setProperties();

    sc = service("FluxSvc",m_fsvc);
    		
    // Open the file and set up the tree:
  
   events = new TTree("Events","Events");
  
	events->Branch("energy",&energy,"energy/D");
	events->Branch("time",&time,"time/D");
	events->Branch("trigger_flag",&trigger_flag,"trigger_flag/I");
	events->Branch("cos_theta",&cos_theta,"cos_theta/D");
	events->Branch("phi",&phi,"phi/D");
	events->Branch("l",&l,"l/D");
	events->Branch("b",&b,"b/D");
	counts=0;
	
   events_recon = new TTree("Events_recon","Events_recon");
	events_recon->Branch("energy",&energy_recon,"energy/D");
	events_recon->Branch("time",&time,"time/D");
	events_recon->Branch("trigger_flag",&trigger_flag,"trigger_flag/I");
	events_recon->Branch("cos_theta",&cos_theta_recon,"cos_theta/D");
	events_recon->Branch("phi",&phi_recon,"phi/D");
	events_recon->Branch("l",&l_recon,"l/D");
	events_recon->Branch("b",&b_recon,"b/D");

    return sc;
    
}

StatusCode TDSReadFluxAlg::execute()
{

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    
    sc = readTDSData();
    return sc;
}

StatusCode TDSReadFluxAlg::readTDSData() {
    // Purpose and Method: Retrieve Monte Carlo data from the TDS.
    MsgStream log(msgSvc(), name());

    static Data record;
    SmartDataPtr<Event::EventHeader> evt(eventSvc(), EventModel::EventHeader);
   
    if(!evt) return StatusCode::SUCCESS;
    //Get the time
    time=(double)evt->time();
    record.time   = time;
    record.energy = 0.0;
    record.theta  = 0.0;
    record.phi    = 0.0;
    
    SmartDataPtr<Event::McParticleCol> particles(eventSvc(), EventModel::MC::McParticleCol);
    if (particles) {
        if(counts % 1000==0) f = new TFile("Events.root","RECREATE"); //or update ?
        log << MSG::INFO << "Retrieved McParticles from the TDS at Time = " 
	<< time <<" s "<< endreq;
	
        Event::McParticleCol::const_iterator p;
	p = particles->begin(); // This is the mother...
	
	energy=(*p)->finalFourMomentum().e();
	dir=(*p)->finalFourMomentum().vect().unit();
	cos_theta=-(*p)->finalFourMomentum().cosTheta(); // note the (-)
	phi=(*p)->finalFourMomentum().phi();
	if(phi < 0) phi += 2*M_PI;
	phi*= deg; //degrees


	trigger_flag=0; // see later...

	Glast2Gal=m_fsvc->transformGlastToGalactic(time);
	GalDirVec=Glast2Gal*dir;
	l = GalDir(GalDirVec).first;
    	b = GalDir(GalDirVec).second;
	//	std::cout<<" ---- Time (sec) "<<time<<std::endl;
	//	std::cout<<"      Energy     = "<<energy<<" (MeV)"<<std::endl;	
	//      std::cout<<"      Direction: Cos(theta) = "<<cos_theta<<std::endl;	
	//	std::cout<<"                        phi = "<<phi<<std::endl;
	//	std::cout<<"-------------------------------------------"<<std::endl;	
 	if(energy>0) events->Fill();
	if(counts % 1000==0) events->Write();
	counts++;
   }
   /////////////////////////////////////////////////////////////////////////////////
    //               		RECONSTRUCTED DATA                                //
    //	 Trigger fag	|   triggered	|    energy>0 	| #track > 0              // 
    // 		0    	|	0	|	0	|	0                 //
    //		1	|	1	|	0	|	0                 //
   //		3	|	1	|	1	|	0                 //
   //		5	|	1	|	0	|	1                 //
   //		7	|	1	|	1	|	1                 //
   /////////////////////////////////////////////////////////////////////////////////
   energy_recon=0.;
   cos_theta_recon=1.;
   phi_recon=0.;
   l_recon=0.;
   b_recon=0.; 

   trigger_flag=0;
   
   SmartDataPtr<Event::CalClusterCol> clusters(eventSvc(),  EventModel::CalRecon::CalClusterCol);
   SmartDataPtr<Event::TkrVertexCol>  tracks(eventSvc(),EventModel::TkrRecon::TkrVertexCol);
   
   if(clusters)	
	{	
	   trigger_flag+=1; // Triggered
//	   cout<<" Cluster Size = "<<clusters->size()<<endl; 
	   energy_recon = clusters->getCluster(0)->getEnergySum();//*1e-3; // convert to GeV
	   log << MSG::INFO << "Reconstructed Energy = " 
	<< energy_recon <<" (MeV) "<< endreq;
	   if (energy_recon>0) trigger_flag+=2; //Triggered + e>0

	if (tracks->size() != 0)
		{
		trigger_flag+=4;

	        const Event::TkrVertex *track = (tracks->front());
        
        	dir = track->getDirection();
        
		cos_theta_recon = -dir.z(); // note the (-)
		phi_recon   = atan2(dir.y(),dir.x());
		if(phi_recon < 0) phi_recon += 2*M_PI;
		phi_recon   *= deg; 
		GalDirVec=Glast2Gal*dir;
		l_recon = GalDir(GalDirVec).first;
    		b_recon = GalDir(GalDirVec).second;	
		log << MSG::INFO << "Reconstructed Direction : Cos Theta = " 
		    << cos_theta_recon <<" phi = "<< phi_recon<<endreq;
		log << MSG::INFO << " Galactic Coordinates : l = " 
		    << l_recon <<" b = "<< b_recon << endreq;
	
	}
	events_recon->Fill(); // The root tree is filled with all the triggered events. 
	// The trigger flag will help to separate them.
   }

   // In the ascii file only the events with energy > 0 and # tracks > 0 will be saved. 
   std::string m_save_opt="energy && track";
   
   bool save=false;
   // This are other options that can be saved
   if(m_save_opt == "energy" && energy_recon>0) 
   	save=true;
   else if(m_save_opt == "track" && tracks->size()>0) 
   	save=true;
   else if(m_save_opt == "triggered" && (clusters || tracks))
   	save=true;
   else if(m_save_opt == "energy || track" && (energy_recon>0 || tracks->size()>0)) 
   	save=true;
   else if(m_save_opt == "energy && track" && (energy_recon>0 &&  tracks->size()>0)) 
   	save=true;
	
   
   if (save)
   {
  	record.energy = energy_recon*1e-3; //Convert to GeV
	record.theta  = acos(cos_theta_recon)*deg;
	record.phi    = phi_recon;
	m_data.push_back(record);
   }
   return  StatusCode::SUCCESS;
}
    

StatusCode TDSReadFluxAlg::finalize()
{
    if(m_save_file.empty()) return StatusCode::SUCCESS;
    std::vector<std::string>::iterator itr;
    for(itr=m_save_file.begin();itr != m_save_file.end();++itr)
    {
      if((*itr)=="root") 
      {
	    f = new TFile("Events.root","RECREATE");
    	    events->Write();
    	    f->Close();
    
    	    f_recon = new TFile("Events_recon.root","RECREATE");
            events_recon->Write();
            f_recon->Close();
      }
      // Write the ascii file
      if((*itr)=="ascii") 
	{
    	long sz = m_data.size();
	if (sz > 0)   
	{
		std::ofstream os("grb_detected.lis");
		os << "nEvents: " << sz << std::endl;
		#ifdef WIN32
		for (std::vector<Data>::iterator it=m_data.begin(); it!=m_data.end(); ++it)
			os << std::setw(12) << std::fixed << it->time << "  " << it->energy << "  " << it->theta << "  " << 
			it->phi << std::endl;
		#else
		for (std::vector<Data>::iterator it=m_data.begin(); it!=m_data.end(); ++it)
			os << std::setw(12) << it->time << "  " << it->energy << "  " << it->theta << "  " << 
			it->phi << std::endl;
		#endif
	}
      }
    }
    return StatusCode::SUCCESS;
}
