#ifndef GRBobsPULSE_HH
#define GRBobsPULSE_HH 1
#include <math.h>
#include "GRBobs/GRBobsConstants.h"
/*! 
 * \class GRBobsPulse
 * \brief GRB phenomenological pulse description
 *
 * The pulse is the elementary structure of a GRB.
 * \author Nicola Omodei       nicola.omodei@pi.infn.it 
 *  
 */


class GRBobsPulse
{
 public:
  GRBobsPulse();
  GRBobsPulse(double peakTime, GRBobsParameters *params);
  void Print();
  inline void SetPeakTime(double t)    {m_peakTime  = t;}
  inline void SetRiseTime(double t)    {m_riseTime   = t;}
  inline void SetDecayTime(double t)   {m_decayTime = t;}
  inline void SetPeakedness(double a)  {m_Peakedness  = a;}
  inline void SetIntensity(double i)   {m_Intensity  = i;}
  
  inline double GetPeakTime()     {return m_peakTime ;}
  inline double GetEndTime()      {return m_end;}
  inline double GetStartTime()    {return m_start;}
  inline double GetDuration()     {return m_duration;}
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
    \\   
    \sigma_r(E)=\sigma_r*\frac{E}{E0}^{-We} \\
    \\
    \sigma_d(E)=\sigma_d*\frac{E}{E0}^{-We} \\
    \end{array}
    \right. 
    \f]
    The variables \f$\sigma_r\f$, \f$\sigma_d\f$ and \f$\nu\f$ are computed by GRBobsParameters::GenerateParameters(), 
    and the constants \em E0,\em We are defined in the namespace ObsCst (ObsCst::E0, ObsCst::We).
    The peaks at different energies are also shifted, so that:
    \f[
    t_p(E0)-t_p(E)= \Delta_t~\sigma_r(1-\frac{E}{E0}^{-We})(log(100.))^{1.0/\nu}
    \f]
    with \f$\Delta_t\f$ corresponding to ObsCst::deltaTPeak.
  */
  double PulseShape(double t ,double e);
 
 private:
  double m_peakTime;
  double m_riseTime;
  double m_decayTime;
  double m_start;
  double m_end;
  double m_duration;
  double m_Intensity;
  double m_Peakedness;
  double m_Epeak;
  double m_LowEnergy;
  double m_HighEnergy;

};

#endif
