// Include files
#include "FluxSvc/IFluxSvc.h"
#include "flux/IFlux.h"
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
//#include "other/GRBTest.cxx"
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
};


static const AlgFactory<GRBTestAlg>  Factory;
const IAlgFactory& GRBTestAlgFactory = Factory;

//------------------------------------------------------------------------------
//
GRBTestAlg::GRBTestAlg(const std::string& name, ISvcLocator* pSvcLocator) :
  Algorithm(name, pSvcLocator){
}

//------------------------------------------------------------------------------
StatusCode GRBTestAlg::initialize() {
  
  
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "initializing..." << endreq;
  // Use the Job options service to set the Algorithm's parameters
  setProperties();
  return StatusCode::SUCCESS;
}


//------------------------------------------------------------------------------
StatusCode GRBTestAlg::execute() 
{
  return StatusCode::SUCCESS;
}


//------------------------------------------------------------------------------
StatusCode GRBTestAlg::finalize() 
{
  return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------------
