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
  m_energy.clear();	
  m_de.clear();  
  m_energy = m_spectobj.getEnergyVector();
  m_de     = m_spectobj.getBinVector();
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
  m_energy.clear();	
  m_de.clear();  
  m_energy = m_spectobj.getEnergyVector();
  m_de     = m_spectobj.getBinVector();
  //////////////////////////////////////////////////
  m_synchrotron = GRBSynchrotron(m_spectobj);
  
  //////////////////////////////////////////////////
}	

/*------------------------------------------------------*/
GRBSim::~GRBSim()
{
  delete m_engine;
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
  m_engine = GRBConstants::GetTheRandomEngine(m_seed);
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
std::vector<double> GRBSim::ComputeFlux(const double time)
{
  // Compute the flux at time t, and return a vector containing the flux in 
  //-------------->       ph/s/MeV/m^2
  //  std::cout<<"Compute the flux @ time ="<<time<<std::endl;
  m_spectobj *= 0.0 ;
  std::vector<double> spectrum(enstep+1 ,0.);
  if (time<0) return spectrum;
  std::vector<GRBShock>::iterator itr=theShocks.begin();
  while(itr != theShocks.end())
    {
      
      m_synchrotron.load((itr),time,m_jetangle,m_distance);
      m_spectobj+=m_synchrotron.getSpectrumObj();
      
      if (cst::flagIC!=0.0)
	{
	  m_icompton = GRBICompton(m_synchrotron.getSpectrumObj());
	  m_icompton.load((itr),time,m_jetangle,m_distance);
	  //m_icompton.setSpectrumObj(m_synchrotron.getSpectrumObj());
	  m_spectobj+=m_icompton.getSpectrumObj();
	}
      
      itr++;
    }
  m_spectobj/=m_area;  // is in ph/s/MeV/m^2
  return m_spectobj.getSpectrumVector();  
}

/*------------------------------------------------------*/
double GRBSim::IFlux(std::vector<double> spctrmVec,double enmin,double enmax)
{ 
  if(spctrmVec.size()==0) return 0.0;
  std::vector<double>::iterator it_spec;
  std::vector<double>::iterator it_ene;
  std::vector<double>::iterator it_de;
  // The SpectrumVec is in ph/s/MeV/m^2  
  // flux is in eV/s/m^2
  double flux=0.0;
  it_ene=m_energy.begin();
  it_spec=spctrmVec.begin();
  for(it_de= m_de.begin();it_de != m_de.end();++it_de)
    {
      if((*it_ene)>=enmin && (*it_ene)<enmax)
	{  
	  flux += (*it_ene)*(*it_spec)*(*it_de)*(1.0e-6);
	  //eV/s/m^2
	}
      it_ene++;
      it_spec++;
    }
  return flux;  
}

/*------------------------------------------------------*/
double GRBSim::IRate(std::vector<double> spctrmVec,double enmin,double enmax)
{
  if(spctrmVec.size()==0) return 0.0;
  std::vector<double>::iterator it_spec;
  std::vector<double>::iterator it_ene;
  std::vector<double>::iterator it_de;
  
  // rate of particle arrivals for enmin <= energy < enmax
  double rate=0.0;
  it_ene=m_energy.begin();
  it_spec=spctrmVec.begin();
  for(it_de= m_de.begin();it_de != m_de.end();++it_de)
    // for (int en=0;en<enstep;en++)
    {
      if((*it_ene)>=enmin && (*it_ene)<enmax)
	{  
	  rate += (*it_spec)*(*it_de)*(1.0e-6);
	  //ph/s/m^2
	}
      it_ene++;
      it_spec++;
    }
  return rate;  
}
/*------------------------------------------------------*/

double GRBSim::DrawPhotonFromSpectrum(std::vector<double> spctrmVec, 
	float u, double enmin)
{
  // With std:iterator
  std::vector<double>::iterator it_spec = spctrmVec.begin();  ;
  std::vector<double>::iterator it_ene;
  std::vector<double> Integral;
  
  if(spctrmVec.size()==0) return 0.0;
  //STEP 1: we need to remove low energy part of the spectrum,
  // in order to avoid drawing photons of no interest to GLAST.
  // minbin: ebergy bin # after which the energy of a photon 
  // will be above emin.

  int minbin=0;
  for(it_ene= m_energy.begin();it_ene != m_energy.end();++it_ene)
	{
	if((*it_ene) >= enmin)
		Integral.push_back((*it_spec)*(*it_ene-*(it_ene-1)));
	else
		minbin++;	
	++it_spec;
	}
  if(u==0) return m_energy[minbin]*1.0e-9;
  //STEP 2:Compute cumulative sum, and normalize to 1.

  for(it_spec=Integral.begin()+1;it_spec!=Integral.end();++it_spec) 
    { 
     *(it_spec-1)-=Integral.front();
      (*it_spec) += *(it_spec-1); 	//Computing cumulative sum
    }

  if( Integral.back() <=0 ) return 0.0;	
  for(it_spec=Integral.begin();it_spec!=Integral.end();++it_spec) 
    { 
      (*it_spec) /= Integral.back(); 	//Normalizing to one
    }
  
  if(Integral.back()!=1.0) return 0.0;	//This should never happen
    
//  std::vector<double>::iterator closer_u = 
//       find_if(Integral.begin(), Integral.end(), bind2nd(greater<double>(), u));
//  --closer_u;

  //STEP 3: Find in the cumulative sum vector, the bin for which
  //the flat random variable u is closest to the value of Integral

  int nabove = Integral.size();
  int nbelow = 0;
  int ibin   = 0;
  while(nabove-nbelow > 1) 
    {
      int middle = (nabove+nbelow)/2;
      if (u == Integral[middle-1]) 
        {
          ibin = middle-1;
          break;
        }
      if (u  < Integral[middle-1]) nabove = middle;
      else                         nbelow = middle;
      ibin = nbelow-1;
    }
  //STEP4: retrurns the centered value of the energy at bin position
  // determined at STEP 3 
    
  double ph = m_energy[ibin+minbin]+
    (m_energy[ibin+minbin+1]-m_energy[ibin+minbin])*
    (Integral[ibin+1] - u)/(Integral[ibin+1] - Integral[ibin]);
    
  return ph*1.0e-9; //returns value in GeV
}

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
