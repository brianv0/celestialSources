#ifndef GRBmanager_H
#define GRBmanager_H

#include <vector>
#include <string>
#include <map>
#include <cmath>
#include "FluxSvc/ISpectrum.h"
#include "facilities/Observer.h"
#include "src/GPS.h"
#include "CLHEP/Random/RandomEngine.h"
#include "GRBSpectrum.h"

class ISpectrum;

/*!\class  GRBmanager
  
  \brief Spectrum class for many GRBs inheriting from GRBSpectrum.
  This class concatenates several GRBSpectrum one after the other 
  for simulating a series of several GRBs.
  
  \author Nicola Omodei       nicola.omodei@pi.infn.it 
  \author Johann Cohen-Tanugi johann.cohen@pi.infn.it
*/
class GRBmanager : public ISpectrum
{

 public:
  /*! This initialize the simulation parseing the parameters.
    
    \param params are set in the xml source library in xml directory.
    They are: 
    - The time of the first burst
    - The time to whait between the next burst

    An example the xml source declaration for this spectrum should appears:
    \verbatim
    <source name=" GRBmanager_Gal">
    <spectrum escale="GeV"> <SpectrumClass name="GRBmanager" params="50 100"/>
    <use_spectrum frame="galaxy"/> 
    </spectrum> </source>
    \endverbatim
  */
  
  GRBmanager(const std::string& params);
  
  virtual  ~GRBmanager(); 
  /*! If a burst is shining it returns the \e GRBSpectrum::flux(time) method 
   */
  double flux(double time)const;
  /* \brief Returns the time interval
   *
   * If a burst is shining it returns the GRBSpectrum::interval(time) method.
   * If not it returns the time to whait for the first photon of the next burst.
   */
  double interval(double time);
  //! \retval 1
  inline double solidAngle() const{return 1.0;}
  //! direction, taken from GRBSim
  inline std::pair<double,double> 
    dir(double energy, HepRandomEngine* engine)
    {return m_GRB->dir(energy, engine);} 
  
  float operator() (float u) const;
  double energySrc(HepRandomEngine*, double time);
  
  std::string title() const {return "GRBmanager";} 
  const char * particleName() const {return "gamma";}
  const char * nameOf() const {return "GRBmanager";}
  
  /*! 
    This method is used to parse the parametyer list
    \param input is the string to parse
    \param index if the position of the parameter in the input list. 
    \retval output is the value of the parameter as float number.
  */  
  float parseParamList(std::string input, int index);  
  
 private:
  
  GRBSpectrum* m_GRB;
  const std::string& m_params;
  float m_FirstTime;
  float m_initialTime;
  float m_timeToWait;
  float m_endTime;
};
#endif
