#include <iterator>
#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <cassert>
#include <string>

#include "GRBSim.h"

/*------------------------------------------------------*/
using namespace cst;
using namespace std;
/*------------------------------------------------------*/

GRBSim::GRBSim(long seed)
{
  std::cout<<"******Staring The GRB Simulation******"<<std::endl;
  m_seed = seed;
  //  std::cout<<"--------INT SEED = "<<m_seed<<std::endl;
  //////////////////////////////////////////////////
  m_spectobj = SpectObj(cst::enmin,cst::enmax,cst::enstep);
  //////////////////////////////////////////////////
  //m_energy.clear();	
  //m_de.clear();  
  //m_energy = m_spectobj.getEnergyVector();
  //m_de     = m_spectobj.getBinVector();
  //////////////////////////////////////////////////
  m_synchrotron = GRBSynchrotron(m_spectobj);
  
  //////////////////////////////////////////////////
}	

GRBSim::GRBSim(const std::string& params)
{
  std::cout<<"******Staring The GRB Simulation******"<<std::endl;
  std::ifstream is (params.c_str());
  is>>m_seed;
  is.close();
  m_seed++;
  std::ofstream os (params.c_str());
  os<<m_seed;
  os.close();
  //  std::cout<<"------PARAMS  SEED = "<<m_seed<<std::endl;
  //////////////////////////////////////////////////
  m_spectobj = SpectObj(cst::enmin,cst::enmax,cst::enstep);
  //////////////////////////////////////////////////
  //m_energy.clear();	
  //m_de.clear();  
  //m_energy = m_spectobj.getEnergyVector();
  //m_de     = m_spectobj.getBinVector();
  //////////////////////////////////////////////////
  m_synchrotron = GRBSynchrotron(m_spectobj);
  
  //////////////////////////////////////////////////
}	

/*------------------------------------------------------*/
GRBSim::~GRBSim()
{
  //  delete m_engine;
  delete myParam;
  std::cout<<"*******Exiting The GRB Simulation ******"<<std::endl;
}
/*------------------------------------------------------*/

void GRBSim::MakeGRB(double time_offset)
{
  //////////////////////////////////////////////////
  m_duration=0.0;
  //////////////////////////////////////////////////
  theShocks.clear();
  try
    {
      myParam=new GRBConstants(); //<- Read the Param
    }
  catch (char * )
    {
      std::cout<< "Failure initializing the GRB constants \n";
      exit(1);
    }
  
  m_enph=myParam->EnergyPh();
  // This initialize the Random Engine in GRBSim...
  //m_engine = GRBConstants::GetTheRandomEngine(m_seed);
  //  GRBengine * m_CentralEngine = GRBengine(myParam);
  GRBengine *m_CentralEngine = new GRBengine(myParam);
  theShocks   = m_CentralEngine->getShocksVector();
  m_duration  = m_CentralEngine->getDuration();
  m_distance  = m_CentralEngine->getDistance();
  m_direction = m_CentralEngine->getDirection();
  m_area      = (4.*M_PI)*pow(m_distance,2.)*1.0e-4; // [m^2]
  m_jetangle  = myParam->JetAngle();
  delete m_CentralEngine;
 
  //////////////////////////////////////////////////
  // All is ok
  ////////////////////////////////////////////////// 
  
  myParam->Print();
  //myParam->Save(cst::savef);

  
  //////////////////////////////////////////////////
  std::cout<<" Dist  of the source   = "<<m_distance<<" cm "<<std::endl;
  std::cout<<" Galactic Direction = [l= "<<m_direction.first<<" ,b = "<<m_direction.second<<"]"<<std::endl;
  std::cout<<" Number of Shocks      = " <<theShocks.size()<< std::endl;
  std::cout<<" Duration of the Burst = "<<m_duration-time_offset<<std::endl;
  //std::cout<<" Total Energy Radiated = "<<0<<std::endl;
  //////////////////////////////////////////////////
  //  std::cout<<m_energy.size()<<" "<<m_de.size()<<std::endl;
  if (m_duration-time_offset<2.)
    {std::cout<<"The burst is Short "<<std::endl;}
  else
    {std::cout<<"The burst is Long "<<std::endl;}
  //ouble dt=(m_duration-time_offset)/nstep;
  
  std::vector<GRBShock>::iterator itr;
  for(itr=theShocks.begin();itr != theShocks.end();++itr)
    {
      (*itr).Write();
    }
}

/*------------------------------------------------------*/
SpectObj GRBSim::ComputeFlux(const double time)
{
  m_spectobj *= 0.0 ;
  if (time<0) return m_spectobj;
  std::vector<GRBShock>::iterator itr=theShocks.begin();
  while(itr != theShocks.end())
    {
      
      m_synchrotron.load((itr),time,m_jetangle,m_distance);
      m_spectobj+=m_synchrotron.getSpectrumObj();
      
      if (cst::flagIC!=0.0)
	{
	  m_icompton = GRBICompton(m_synchrotron.getSpectrumObj());
	  m_icompton.load((itr),time,m_jetangle,m_distance);
	  m_spectobj+=m_icompton.getSpectrumObj();
	}
      
      itr++;
    }
  m_spectobj/=m_area;  // is in ph/s/MeV/m^2
  return m_spectobj;  
}

/*------------------------------------------------------*/
long GRBSim::parseParamList(std::string input, int index)
{
  std::vector<long> output;
  unsigned int i=0;
  for(;!input.empty() && i!=std::string::npos;){
    float f = ::atoi( input.c_str() );
    output.push_back(f);
    i=input.find_first_of(",");
    input= input.substr(i+1);
  } 
  return output[index];
}
