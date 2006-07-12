#ifndef GRBobsPULSE_HH
#define GRBobsPULSE_HH 1
#include <math.h>
/*! 
 * \class GRBobsPulse
 * \brief GRB phenomenological pulse description
 *
 * The pulse is the elementary structure of a GRB
 * \author Nicola Omodei       nicola.omodei@pi.infn.it 
 *  
 */


class GRBobsPulse
{
 public:
  GRBobsPulse();
  GRBobsPulse(double, double, double, double, double, double, double, double);
  void Print();
  inline void SetPeakTime(double t)    {m_peakTime  = t;}
  inline void SetRiseTime(double t)    {m_riseTime   = t;}
  inline void SetDecayTime(double t)   {m_decayTime = t;}
  inline void SetPeakedness(double a)  {m_Peakedness  = a;}
  inline void SetIntensity(double i)   {m_Intensity  = i;}
  
  inline double GetPeakTime()     {return m_peakTime ;}
  inline double GetEndTime()      {return m_peakTime + log(100.0) * m_decayTime;}
  inline double GetStartTime()    {return m_peakTime - log(100.0) * m_riseTime;}

  /*!  
    The pulse shape is composed by a temporal profile of equation:
    \f[
    \left\{
    \begin{array}{ll}
    I \exp[-(|t-t_{peak} |/\sigma_r)^\nu],  \leq t_{peak}\\
    \\
    I \exp[-(|t-t_{peak} |/\sigma_d)^\nu], t>t_{peak}\\
    \end{array}\right.
    \f]
    where \f$ I \f$ is the intensity of the pulse,\f$t_{peak}\f$ is the peak time,\f$ \sigma_r\f$ and \f$ \sigma_d\f$ are the rise time and the decay time of the pulse,\f$ \nu \f$ is the peakedness.
    The spectral shape is the Band function:
    \f[
    N(E) = N_0
    \left\{
    \begin{array}{l l}
    (E)^{\alpha} \exp (-{E\over E_0}), & for E < (\alpha-\beta)E_0 \\
    \\
    (\alpha-\beta) E_0^{(\alpha-\beta)} (E)^{\beta} \exp(\beta-\alpha), & for E
    > (\alpha-\beta)E_0,
    \end{array}
    \right. 
    \f]
    
    Other relation used for building the model are:
    \f[
    \left\{
    \begin{array}{l}
    W(E) = W_0 E^{-0.4} \\
    \\   
    \sigma_r=0.33\sigma_d^{0.83}\\
    \\
    \sigma_d = 0.75~0.69^{1/\nu}~W\\
    \end{array}
    \right. 
    \f]
    The first equation expresses the pulse width dependence on the energy. The second equation correlates the rise time with the decay time, and the third correlates the width with the decay time.    
   */
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
