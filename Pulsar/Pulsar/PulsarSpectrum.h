#ifndef PulsarSpectrum_H
#define PulsarSpectrum_H

#include <vector>
#include <string>
#include <map>
#include <cmath>
#include "flux/ISpectrum.h"
#include "flux/EventSource.h"
#include "facilities/Util.h"
#include "astro/JulianDate.h" 
#include "PulsarConstants.h"
#include "PulsarSim.h"
#include "SpectObj/SpectObj.h"

/*! 
* \class PulsarSpectrum
* \brief Class that starts the Pulsar simulation according to the parameters specified in the XML file. 
*  
* \author Nicola Omodei        nicola.omodei@pi.infn.it 
* \author Massimiliano Razzano massimiliano.razzano@pi.infn.it
*
* It takes the parameters of the model and the name of the pulsar from the XML file. Then looks in the DataList file for the* name of the pulsar and extract the speficic parameters of the pulsar (period, flux, ephemerides, etc.). 
* The DataFile is placed in <i>data</i> directory. This class, derived from ISpectrum, takes into account 
* the decorretions for the ephemerides, as described in the method <b>Interval</b>.
*/

class PulsarSpectrum : public ISpectrum
{
  
 public:

  //! Constructor of PulsarSpectrum  
  PulsarSpectrum(const std::string& params);
  
  //! Constructor of PulsarSpectrum  
  virtual  ~PulsarSpectrum();
   
  //! Return the flux of the Pulsar at time t
  double flux(double time)const;

  //! Returns the time interval to the next photon
  double interval(double time);

  //! Returns the number of turns made at a specified time
  double getTurns(double time);

  //Retrieve the nextTime from the number of turns to be completed
  double retrieveNextTimeTilde( double tTilde, double totalTurns, double err );
  
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

  //! Parse parameters from XML file
  std::string parseParamList(std::string input, unsigned int index);  
      
 private:
  
  PulsarSim *m_Pulsar;
  SpectObj  *m_spectrum;
  
  astro::JulianDate *JDStartMission; 
  astro::JulianDate *JDStartSimulation; 
  astro::JulianDate *m_JDCurrent; 

  const std::string& m_params; 
  
  std::string m_PSRname;
  double m_RA, m_dec, m_l, m_b;  
  double m_period, m_pdot, m_p2dot, m_t0, m_phi0, m_f0, m_f1, m_f2;
  int m_numpeaks;
  int m_model;
  double m_flux, m_enphmin, m_enphmax;
  int m_seed;
  
};
#endif