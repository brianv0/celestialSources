#include <iostream>
#include <algorithm>

#include "GRBobs/GRBobsengine.h"
#include "GRBobs/GRBobsPulse.h"

#define DEBUG 0

using namespace ObsCst;

GRBobsengine::GRBobsengine(GRBobsParameters *params)
  : m_params(params)
{
  m_dir = m_params->GetGalDir(); 
}

double GRBobsengine::generatePulses(std::vector<GRBobsPulse*> &thePulses, double duration)
{
  //////////////////////////////////////////////////
  double tau, pt1,pt,rt, dt, ph, nu, ep, a, b,endTime,BurstEndTime;

  pt=0.0;
  
  GRBobsPulse *aPulse;
  endTime=0.0;
  BurstEndTime=0.0;
  int npulses=0;
  
  while(endTime<duration || npulses==0)
    {
      m_params->GenerateParameters();
      tau = m_params->GetPulseSeparation();
      rt  = m_params->GetRiseTime();
      dt  = m_params->GetDecayTime();
      ph  = m_params->GetPulseHeight();
      nu  = m_params->GetPeakedness();
      ep  = m_params->GetEpeak();
      a   = m_params->GetLowEnergy();
      b   = m_params->GetHighEnergy();
      if (npulses==0) 
	pt  = pow(log(100.0),1./nu) * rt; //this sets the tstart =0
      else 
	pt=pt1+tau; 
      
      aPulse = new GRBobsPulse(pt,rt,dt,ph,nu,ep,a,b);
      
      if(DEBUG) 
	{
	  m_params->PrintParameters();
	  aPulse->Print();
	}

      endTime = aPulse->GetEndTime();
      
      if(endTime <= duration) 
	{
	  thePulses.push_back(aPulse);
	  pt1 = pt;
	  npulses++;
	  BurstEndTime = TMath::Max(BurstEndTime,endTime);
      	}
      else
	{
	  delete aPulse;
	}
    }
  
  return endTime;
}

std::vector<GRBobsPulse*> GRBobsengine::CreatePulsesVector()  
{
  // get goal duration value;
  double duration = m_params->GetDuration();

  double burstEndTime;
  std::vector<GRBobsPulse*> thePulses;
  bool done = false;

  do {
    burstEndTime = generatePulses(thePulses, duration);
    if (burstEndTime < 0.99 * duration) 
      {    // return resources and start over
	for (unsigned i = 0; i < thePulses.size(); i++) delete thePulses[i];
	thePulses.erase(thePulses.begin(), thePulses.end());
      }
    else done = true;
  }  while (!done) ;
  return thePulses;
}

//////////////////////////////////////////////////

/*
  double GRBobsengine::getDistance()
  {
  double m_redshift = cst::red;
  double qo=(1.0+3.0*cst::wzel)/2.0;
  return ((cst::c/(cst::Hubble*1.0e+5)/pow(qo,2.0))*
  (m_redshift*qo+(qo-1.0)*(-1.0+sqrt(2.0*qo*m_redshift+1.0)))*cst::mpc2cm);
  }
*/


