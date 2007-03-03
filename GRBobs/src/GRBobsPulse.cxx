#include "GRBobs/GRBobsPulse.h"
#include "GRBobs/GRBobsConstants.h"
#include <iostream>
#include <cmath>

#define DEBUG 0

using std::fabs; using std::pow;

GRBobsPulse::GRBobsPulse(){;}

GRBobsPulse::GRBobsPulse(double peakTime,
			 double riseTime,
			 double decayTime,
			 double Intensity,
			 double Peakedness,
			 double Epeak,
			 double LowEnergy,
			 double HighEnergy)
{  
  m_peakTime     = peakTime;
  m_riseTime     = riseTime;
  m_decayTime    = decayTime;
  m_Intensity    = Intensity;
  m_Peakedness   = Peakedness;
  m_start        = peakTime - riseTime  * pow(log(100.0),1.0/Peakedness);
  m_end          = peakTime + decayTime * pow(log(100.0),1.0/Peakedness);
  m_Epeak        = Epeak;
  m_LowEnergy    = LowEnergy;
  m_HighEnergy   = HighEnergy;
  m_duration     = m_end-m_start;
  if(DEBUG)  Print();
  
}

void GRBobsPulse::Print()
{
  std::cout<<"Pulse: start= "<<m_start
	   <<" peak= "<<m_peakTime
	   <<" end= "<<m_end
	   <<" rise t= "<<m_riseTime
	   <<" Decay t= "<<m_decayTime
	   <<" Dur.= "<<m_duration
	   <<" Intensity= "<<m_Intensity
	   <<" Peakedness= "<<m_Peakedness
	   <<" Peak E= "<<m_Epeak
	   <<" Low idx= "<<m_LowEnergy-ObsCst::We
	   <<" High idx= "<<m_HighEnergy-ObsCst::We
	   <<std::endl;
}

double GRBobsPulse::PulseShape(double t, double e)
{
  //  double tp = m_ts + m_tp;
  double rt = m_riseTime  * pow(e/ObsCst::E0,-ObsCst::We);
  double dt = m_decayTime * pow(e/ObsCst::E0,-ObsCst::We);
  double deltaTP = ObsCst::deltaTPeak * (m_riseTime - rt) * pow(log(100.),1.0/m_Peakedness);
  
  double pt = m_peakTime - deltaTP;
  double pulse=0;
  if (t<pt)
    {
      pulse = exp(-pow(fabs(t-pt)/rt,m_Peakedness));
    }
  else
    {
      pulse = exp(-pow(fabs(t-pt)/dt,m_Peakedness));
    }
  
  //////////////////////////////////////////////////
  double a  = m_LowEnergy;
  double b  = m_HighEnergy;
  
  double Eb = m_Epeak;
  double Ep = (2.+a)/(a-b)*Eb;
  double C  = pow(Eb,a-b)*exp(b-a);
  
  double bandf;
  if(e < Eb) 
    bandf = pow(e,a) * exp(-e*(2.+a)/Ep);
  else
    bandf= C * pow(e,b); // ph cm^(-2) s^(-1) keV^(-1)
  
  return m_Intensity * bandf * pulse;
}



