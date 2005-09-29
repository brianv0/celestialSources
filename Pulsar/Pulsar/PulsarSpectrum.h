/////////////////////////////////////////////////
// File PulsarSpectrum.h
// Header file for PulsarSpectrum class
//////////////////////////////////////////////////

#ifndef PulsarSpectrum_H
#define PulsarSpectrum_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>
#include <map>
#include <cmath>
#include "PulsarConstants.h"
#include "PulsarSim.h"
#include "SpectObj/SpectObj.h"
#include "flux/ISpectrum.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "flux/EventSource.h"
#include "facilities/Util.h"
#include "astro/JulianDate.h"
#include "astro/EarthOrbit.h"
#include "astro/SolarSystem.h"
#include "astro/GPS.h"

/*! 
* \class PulsarSpectrum
* \brief Class that generates a Pulsar sources according to the parameters specified in the XML file (model parameters) and in 
*        the PulsarDataList.txt located in the /data directory. 
*  
* \author Nicola Omodei        nicola.omodei@pi.infn.it 
* \author Massimiliano Razzano massimiliano.razzano@pi.infn.it
*
* PulsarSpectrum, derived from ISpectrum, takes the model parameters from the XML file, where are also located the pulsar name, his position in the sky, 
* and the energy range of the extracted photons. Then it looks in the PulsarDataList.txt file ( in the <i>/data</i> directory)
* for the name of the pulsar and then extracts the specific parameters of the pulsar (period, flux, ephemerides, etc.) related 
* to that pulsar. 
* Then it computes all the needed decorretions for timing, in particular the period changes and the barycentric decorrections.
*/

class PulsarSpectrum : public ISpectrum
{
  
 public:

  //! Constructor of PulsarSpectrum  
  PulsarSpectrum(const std::string& params);
  
  //! Destructor of PulsarSpectrum  
  virtual  ~PulsarSpectrum();
   
  //! Return the flux of the Pulsar at time t
  double flux(double time)const;

  //! Returns the time interval to the next extracted photon, according to the flux, in a frame where pdot=0 and no barycentric decorrections
  double interval(double time);

  //! Returns the number of turns made by the pulsar at a specified time, referred to an inital epoch t0
  double getTurns(double time);

  //!Retrieve the nextTime from the number of turns to be completed
  double retrieveNextTimeTilde( double tTilde, double totalTurns, double err );

  //!Apply the barycentric corrections and returns arrival time in TDB
  double getBaryCorr( double tdbInput ); 

  //! Get the decorrected time in TDB starting from a TT corrected time
  double getDecorrectedTime( double CorrectedTime);
  
  //! Get the pulsar ephemerides and data from the DataList
  int getPulsarFromDataList(std::string sourceFileName);

  //! Save an output txt file compatible with D4 file
  int saveDbTxtFile();

  //! direction, taken from PulsarSim
  inline std::pair<double,double>
    dir(double energy) 
    {
      return m_GalDir;
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
  
  astro::EarthOrbit *m_earthOrbit;
  astro::SolarSystem m_solSys;
 
  astro::SkyDir m_PulsarDir;
  Hep3Vector m_PulsarVectDir;

  std::pair<double,double> m_GalDir;
  
  const std::string& m_params; 
  
  std::string m_PSRname;

  double m_RA, m_dec, m_l, m_b;  
  double m_period, m_pdot, m_p2dot, m_t0, m_t0Init, m_t0End, m_phi0, m_f0, m_f1, m_f2;
  std::string m_ephemType;
  std::vector<double> m_periodVect, m_pdotVect, m_p2dotVect, m_f0Vect, m_f1Vect, m_f2Vect, m_phi0Vect, m_t0Vect, m_t0InitVect, m_t0EndVect, m_txbaryVect;
  
  int m_numpeaks;
  int m_model;
  double m_flux, m_enphmin, m_enphmax;
  int m_seed;


};
#endif
