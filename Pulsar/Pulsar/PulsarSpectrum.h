#ifndef PulsarSpectrum_H
#define PulsarSpectrum_H
#include "PulsarConstants.h"

#include <vector>
#include <string>
#include <map>
#include <cmath>
#include "flux/ISpectrum.h"
#include "PulsarSim.h"
#include "SpectObj/SpectObj.h"

#include "facilities/Util.h"

//class ISpectrum;

class PulsarSpectrum : public ISpectrum
{
  
 public:
  /*! This initializes the simulation parsing the parameters.

  */
  
  PulsarSpectrum(const std::string& params);
  
  virtual  ~PulsarSpectrum();
   
  /*! If a burst is shining it returns the PulsarSpectrum::flux method 
   */
  double flux(double time)const;
  /*! \brief Returns the time interval
   *
   * If a burst is shining it returns the PulsarSpectrum::interval method.
   * If not it returns the time to whait for the first photon of the next burst.
   */
  double interval(double time);
  
  //! direction, taken from PulsarSim
  inline std::pair<double,double>
    dir(double energy) 
    {
      return std::make_pair(1,1);
    } 
  //! calls PulsarSpectrum::energySrc
  double energy(double time);
  
  std::string title() const {return "PulsarSpectrum";} 
  const char * particleName() const {return "gamma";}
  const char * nameOf() const {return "PulsarSpectrum";}
  
  /*! 
    This method parses the parameter list
    \param input is the string to parse
    \param index if the position of the parameter in the input list. 
    \retval output is the value of the parameter as float number.
  */  
  double parseParamList(std::string input, int index);  
  
 private:
  
  PulsarSim   *m_Pulsar;
  SpectObj *m_spectrum;


  const std::string& m_params;

  double m_period,m_fluence;
};
#endif