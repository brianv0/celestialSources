// Include files
#include "FluxSvc/IFluxSvc.h"
#include "FluxSvc/IFlux.h"
#include <fstream>
// GlastEvent for creating the McEvent stuff
#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/MCEvent.h"
#include "Event/MonteCarlo/McParticle.h"
#include "Event/TopLevel/EventModel.h"

// Gaudi system includes
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IParticlePropertySvc.h"
#include "GaudiKernel/SmartRefVector.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include <list>
#include <string>
#include <vector>
//#include "GaudiKernel/ParticleProperty.h"


// GRB includes
//#include "../GRB/GRBSpectrum.h"
#include "other/GRBTest.cxx"
//include files for ROOT...
#include "TTree.h"
#include "TBranch.h"
#include "TObjArray.h"
#include "TFile.h"


using namespace std;

/*!\class GRBTestAlg
  \brief Test the GRB Gaudi algorithm.
  
  It contains the structure of a general algorithm of Gaudi.
  All the option available are declered in the joboptions file.
  To test the algorithm, run: 
  \verbatim
  test_GRB.exe ../src/test/GRBtestAlgOptions.txt
  \endverbatim
  
*/
class GRBTestAlg : public Algorithm {
  
public:
  /*! \brief Constructor of this form must be provided

    In the constructor the algorithm declares its properties
   */ 
  GRBTestAlg(const std::string& name, ISvcLocator* pSvcLocator);   
  
  /*! \brief Initializes the algorithm.
    
    It uses the Job options service to set the Algorithm's parameters, 
    and points to flux service 
  */

  StatusCode initialize();
  
  //! Executes the algorithm, calling GRBTest.
  StatusCode execute();
  
  //!Finalizes the algorithm.
  StatusCode finalize();
  
private:
  //! pointer to the flux Service 
  IFluxSvc* m_fsvc;
    
  int m_loop;
  double m_time;
  int m_events;

  std::string m_source_name;
  
  std::vector<std::string> m_background_name;
  std::vector<std::string> m_save_file;
  double m_observation_time;
  
};


static const AlgFactory<GRBTestAlg>  Factory;
const IAlgFactory& GRBTestAlgFactory = Factory;

//------------------------------------------------------------------------------
//
GRBTestAlg::GRBTestAlg(const std::string& name, ISvcLocator* pSvcLocator) :
  Algorithm(name, pSvcLocator),m_time(10.0),m_events(100000),
  m_source_name("GRBSpectrum"){
  
  declareProperty("source_name", m_source_name);
  declareProperty("background_name",  m_background_name);
  declareProperty("observation_time", m_time);
  declareProperty("EvtMax",   m_events);
  declareProperty("savefile", m_save_file);

}

//------------------------------------------------------------------------------
StatusCode GRBTestAlg::initialize() {
  
  
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "initializing..." << endreq;
  
  // Use the Job options service to set the Algorithm's parameters
  setProperties();
  
  if( m_save_file.empty() )
    {
      m_save_file.push_back("no");
    } 
  
  /*
    if( m_background_name.empty() )
    {
    m_background_name.push_back(" ");
    }
  */
  // get the service
  StatusCode sc = service("FluxSvc", m_fsvc);
  m_loop=0;
  std::ofstream ofs("GRBtemp.txt");
  ofs<<0<<std::endl;
  ofs.close();
  return sc;
}


//------------------------------------------------------------------------------
StatusCode GRBTestAlg::execute() {
  m_loop++;
  cout<<m_loop<<endl;
  StatusCode  sc = StatusCode::SUCCESS;
  MsgStream   log( msgSvc(), name() );    
  
  std::vector<char*> arguments;
  arguments.push_back(" ");
  
  char arg1[50];
  char arg2[50];
  char arg3[50];

  sprintf(arg1,"%f",m_time);
  sprintf(arg2,"%d",m_events);
  sprintf(arg3,"%s%d",m_source_name.c_str(),m_loop);
  sprintf(arg3,"%d",m_loop);
  
  arguments.push_back("-time");
  arguments.push_back(arg1);
  arguments.push_back("-events");
  arguments.push_back(arg2);
  
  std::vector<std::string>::iterator itr;
    
  for(itr=m_save_file.begin();itr != m_save_file.end();++itr)
    {
      if((*itr)=="root") 
	{
	  arguments.push_back("-root");
	  arguments.push_back(arg3);
	  std::cout<<" Saving "<<(arg3)<<" in a root file..."<<std::endl;
	}
      if((*itr)=="ascii") 
	{
	  arguments.push_back("-ascii");
	  arguments.push_back(arg3);
	  std::cout<<" Saving the events in an ascii file..."<<std::endl;
	}
    }
  
  arguments.push_back(const_cast<char *>(m_source_name.c_str()));
  std::cout<<"Primary Source Name = "<<m_source_name.c_str()<<std::endl;
  if( !m_background_name.empty() ){
    for(itr=m_background_name.begin();itr !=m_background_name.end();++itr)
      {
	arguments.push_back(const_cast<char *>((*itr).c_str()));
	std::cout<<" Background Source Name = "<<const_cast<char *>((*itr).c_str())<<std::endl;
      }
  }
  
  
  GRBTest* m_GRBTest = new GRBTest();
  
  m_GRBTest->setService(m_fsvc);
  //m_GRBTest->help();
  //m_GRBTest->listSpectra();
  m_GRBTest->Start(arguments);
  
  return sc;
}


//------------------------------------------------------------------------------
StatusCode GRBTestAlg::finalize() {
    
    return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------------
