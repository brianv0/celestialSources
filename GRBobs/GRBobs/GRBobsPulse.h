#ifndef GRBobsPULSE_HH
#define GRBobsPULSE_HH 1
#include <math.h>
class GRBobsPulse
{
 public:
  GRBobsPulse();
  GRBobsPulse(double, double, double, double, double, double, double, double);
  
  inline void SetPeakTime(double t)    {m_peakTime  = t;}
  inline void SetRiseTime(double t)    {m_riseTime   = t;}
  inline void SetDecayTime(double t)   {m_decayTime = t;}
  inline void SetPeakedness(double a)  {m_Peakedness  = a;}
  inline void SetIntensity(double i)   {m_Intensity  = i;}
  
  inline double GetPeakTime()     {return m_peakTime ;}
  inline double GetEndTime()      {return m_peakTime + log(100.0) * m_decayTime;}
  inline double GetStartTime()    {return m_peakTime - log(100.0) * m_riseTime;}
  double PulseShape(double t ,double e);
 
 private:
  double m_peakTime;
  double m_riseTime;
  double m_decayTime;
  double m_Intensity;
  double m_Peakedness;
  double m_Epeak;
  double m_LowEnergy;
  double m_HighEnergy;

};

#endif
