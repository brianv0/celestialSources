/*!
  \class  GRBmanager
  
  \brief Spectrum class for many GRBs 
 
  This class concatenates several GRB one after the other 
  for simulating a series of several GRBs.
  This class fill the interface provided by ISpectrum class.

  \author Nicola Omodei       nicola.omodei@pi.infn.it 
  \author Johann Cohen-Tanugi johann.cohen@pi.infn.it
*/

#ifndef GRBmanager_H
#define GRBmanager_H
#include "GRBConstants.h"

#include "TString.h"

#include <vector>
#include <string>
#include <cmath>
#include "flux/ISpectrum.h"
#include "flux/EventSource.h"
#include "GRBSim.h"
#include "SpectObj/SpectObj.h"

#include "facilities/Util.h"

class GRBmanager : public ISpectrum
{
  
 public:
  /*! This initializes the simulation parsing the parameters.
    
    \param params are set in the xml source library in xml directory.
    They are: 
    - The time of the first burst
    - The average time to wait before the next burst
    - The Minimum photon energy sampled.
    
    An example the xml source declaration for this spectrum should appears:
    \verbatim
    <source name="GRB">
      <spectrum escale="MeV"> 
        <SpectrumClass name="GRBmanager" params="1000,86400,30.0"/>
        <use_spectrum frame="galaxy"/> 
      </spectrum> 
    </source>
    \endverbatim
    Rapresenting the first burst at 1000 seconds. The intervals fot the other bursts will be sampled from an esponential function with mean 86400 seconds.
    Photons greather than 30 MeV will be sampled and passed to the MonteCarlo (via FluxSvc)
  */
  
  GRBmanager(const std::string& params);
  
  virtual  ~GRBmanager();
  /// It initialize the simulation creating a "burs" object (GRBSim) and a flux object (SpectrObj).
  void GenerateGRB();
  /*! 
    This method parses the parameter list
    \param input is the string to parse
    \param index if the position of the parameter in the input list. 
    \retval output is the value of the parameter as float number.
  */  
  double parseParamList(std::string input, int index);  
  /*! 
    If a burst is ON (shining) it returns the flux of the busrt. It calls the SpectrObj::flux method for computing the flux.
  */
  double flux(double time)const;
  /*! Returns the time interval.
   *
   * If a burst is ON it returns the SpectrObj::interval method for computing the Interval method.
   * If not it returns the time to whait for the first photon of the next burst.
   */
  double interval(double time);

  /// Returns the energy of the photon in MeV (it calls calls Spectrum::energy)
  double energy(double time);
  
  //! direction, taken from GRBSim
  inline std::pair<double,double> dir(double energy) 
    {
      return m_GalDir;
    } 

  std::string title() const {return "GRBmanager";} 
  const char * particleName() const {return "gamma";}
  const char * nameOf() const {return "GRBmanager";}
  /// This method uses astro::EarthOrbit and astro::JulianDate for concatenating the string containing the GRB name, using the default GRB notation.
  /// GRB+YEAR+MONTH+DAT+FRAC_OF_DAY.
  TString GetGRBname(double time);
  
 private:
  double m_Rest;
  double m_Frac;
  double m_ra;
  double m_dec;
  double m_l;
  double m_b;
  double m_theta;
  double m_phi;
  
  GRBSim   *m_GRB;
  SpectObj *m_spectrum;
  Parameters *m_par;
  
  const std::string& m_params;
  std::string paramFile;

  double m_timeToWait;
  double m_startTime;
  double m_fluence;
  UInt_t m_GRBnumber;
  double m_enph;
  double m_endTime;

  double m_nextBurst;
  int    m_Nbursts;
  std::pair<double,double> m_GalDir;
  bool m_GenerateGBMOutputs;
};
#endif
