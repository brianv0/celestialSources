#include "GRBobs/GRBobsPulse.h"
#include <iostream>

#define DEBUG 1  

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
  m_Epeak        = Epeak;
  m_LowEnergy    = LowEnergy;
  m_HighEnergy   = HighEnergy;
  if(DEBUG)  Print();
  
}

void GRBobsPulse::Print()
{
  std::cout<<"Pulse: start= "<<m_peakTime
	   <<" rise t= "<<m_riseTime
	   <<" Decay t= "<<m_decayTime
	   <<" Intensity= "<<m_Intensity
	   <<" Peakedness= "<<m_Peakedness
	   <<" Peak Energy= "<<m_Epeak
	   <<" Low Energy idx= "<<m_LowEnergy-0.4
	   <<" High Energy idx= "<<m_HighEnergy-0.4
	   <<std::endl;
}

double GRBobsPulse::PulseShape(double t, double e)
{
  //  double tp = m_ts + m_tp;
  double rt = m_riseTime  * pow(e/20.0,-0.4);
  double dt = m_decayTime * pow(e/20.0,-0.4);
  double pulse=0;
  if (t<m_peakTime)
    {
      pulse = exp(-pow(fabs(t-m_peakTime)/rt,m_Peakedness));
    }
  else
    {
      pulse = exp(-pow(fabs(t-m_peakTime)/dt,m_Peakedness));
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



