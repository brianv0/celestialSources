#include "GRBobs/GRBobsPulse.h"
#include <iostream>

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
  m_Peakedness   = 1.0;//Peakedness;
  m_Epeak        = Epeak;
  m_LowEnergy    = LowEnergy;
  m_HighEnergy   = HighEnergy;
  std::cout<<GetStartTime()<<" "<<GetPeakTime()<<" "<<GetEndTime()<<std::endl;
}

double GRBobsPulse::PulseShape(double t, double e)
{
  //  double tp = m_ts + m_tp;
  double rt = m_riseTime  * pow(e,-0.4);
  double dt = m_decayTime * pow(e,-0.4);
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
  
  //double Ep = m_Epeak;
  //  double Eb = (a-b)/(2. + a)*Ep;
  double C  = pow(Eb,a-b)*exp(b-a);
  
  double bandf;
  if(e < Eb) 
    bandf = pow(e,a) * exp(-e*(2.+a)/Ep);
  else
    bandf= C * pow(e,b); // ph cm^(-2) s^(-1) keV^(-1)
  
  /*  if(fabs(t-m_peakTime)<0.016 && e < 10.0)
      std::cout<<m_peakTime<<" "<<dt<<" "<<rt<<" "<<m_Peakedness<<" "<<pulse<<" "<<a<<" "<<b<<" "<<Ep<<" "<<Eb<<std::endl;
  */
  return m_Intensity * bandf * pulse;
}



