#include <fstream>
#include <vector>

#include "GaudiKernel/Algorithm.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"

#include "../LatGRBAlert/BackgroundMixer.cxx"
#include "../LatGRBAlert/LatGRBAlert.cxx"
#include "../LatGRBAlert/ClusterData.cxx"


class LatGRBAlertAlg : public Algorithm
{
  public:
	  // Constructor of this form must be provided
	  LatGRBAlertAlg(const std::string& name, ISvcLocator* pSvcLocator);

      // Three mandatory member functions of any Gaudi algorithm
      StatusCode initialize();
      StatusCode execute();
      StatusCode finalize()   { return StatusCode::SUCCESS; }

 private:
	 std::string  m_grbFile;
	 std::string  m_backgroundFile;
	 std::string  m_mixedFile;
	 double       m_grbOffsetTime;
	 long         m_nbckoff;
	 bool         m_mix;
};



// Static Factory declaration
static const AlgFactory<LatGRBAlertAlg>  Factory;
const IAlgFactory& LatGRBAlertAlgFactory = Factory;

LatGRBAlertAlg::LatGRBAlertAlg(const std::string& name, ISvcLocator* pSvcLocator) :
  Algorithm(name, pSvcLocator)
{
	declareProperty("grbFile", m_grbFile); 
	declareProperty("backgroundFile",  m_backgroundFile);
    declareProperty("mixedFile", m_mixedFile);
    declareProperty("grbOffsetTime", m_grbOffsetTime);
    declareProperty("nbckoff", m_nbckoff);
    declareProperty("mix", m_mix);
}



StatusCode LatGRBAlertAlg::initialize() 
{
	MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initializing..." << endreq;
  
    // Use the Job options service to set the Algorithm's parameters
    setProperties();

	return StatusCode::SUCCESS;
}



StatusCode LatGRBAlertAlg::execute()
{
    MsgStream         log( msgSvc(), name() );

	BackgroundMixer *backgroundMix;

	if (m_mix)
	{
		try
		{
			backgroundMix = new BackgroundMixer(m_grbFile, m_backgroundFile, m_mixedFile, 0);
		}

		catch (...)
		{
			log << MSG::ERROR << "Could not allocate memory for BackgroundMixer object." << endreq;
			return StatusCode::FAILURE;
		}
	}

	else
	{
		try
		{
			backgroundMix = new BackgroundMixer(m_mixedFile);
		}

		catch (...)
		{
			log << MSG::ERROR << "Could not allocate memory for BackgroundMixer object." << endreq;
			return StatusCode::FAILURE;
		}
	}

	LatGRBAlert latGRBAlert(backgroundMix->photonData(), m_nbckoff, backgroundMix->ngrb(), backgroundMix->nbck());
	delete backgroundMix;

	std::cout << "like-size: " << latGRBAlert.jointLike().size() << std::endl;

    return StatusCode::SUCCESS;
}