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
#include <stdexcept>
#include "PulsarConstants.h"
#include "PulsarSim.h"
#include "SpectObj/SpectObj.h"
#include "flux/Spectrum.h"
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
* Then it computes all the needed decorretions for timing, in particular:
*  - Period changes according to some ephemerides;
*  - Barycentric corrections:
*  - Timing Noise;
*  - Binary demodulation
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

  //! Returns the time interval to the next extracted photon, according to the flux, in a frame where pdot=0 and no barycentric decorrections and also no binary demodulation
  double interval(double time);

  //! Returns the number of turns made by the pulsar at a specified time, referred to an inital epoch t0
  double getTurns(double time);

  //!Retrieve the nextTime from the number of turns to be completed
  double retrieveNextTimeTilde( double tTilde, double totalTurns, double err );

  //!Apply the barycentric corrections and returns arrival time in TDB
  double getBaryCorr(double tdbInput, int LogCorrFlag);
 
  //! Compute binary demodulation in iterative way
  double getIterativeDemodulatedTime(double tInput, int LogFlag);

  //! Compute binary demodulation in a single step
  double getBinaryDemodulation( double tInput, int LogDemodFlag);

  //! Get the decorrected time in TDB starting from a TT corrected time (inverse of getBaryCorr)
  double getDecorrectedTime( double CorrectedTime);

  //! Get the binary modulated time (inverse of getBinaryDemodulation)
  double getBinaryDemodulationInverse( double CorrectedTime);
  
  //! Get the pulsar ephemerides and data from the DataList
  int getPulsarFromDataList(std::string sourceFileName);

  //! Get the binary pulsar orbital data from the BinDataList
  int getOrbitalDataFromBinDataList(std::string sourceBinFileName);

  //! Save an output txt file with pulsar ephemerides compatible with D4 file
  int saveDbTxtFile();

  //! Save an output txt file with pulsar orbital data compatible with D4 file
  int saveBinDbTxtFile();

  //! Get the eccentric anomaly at time t (only for binaries)
  double GetEccentricAnomaly(double mytime);

  //! calls PulsarSpectrum::energySrc
  double energy(double time);

  //! direction, taken from PulsarSim
  inline std::pair<double,double>
    
    dir(double energy) 
    {
      return m_GalDir;
    } 

  std::string title() const {return "PulsarSpectrum";} 
  const char * particleName() const {return "gamma";}
  const char * nameOf() const {return "PulsarSpectrum";}

  //! Parse parameters from XML file
  std::string parseParamList(std::string input, unsigned int index);  
      
 private:
  
  PulsarSim *m_Pulsar; 
  SpectObj  *m_spectrum;
  TRandom *m_PSpectrumRandom;
  
  //! Variables related to Solar System and Astro.
  astro::EarthOrbit *m_earthOrbit;
  astro::SolarSystem m_solSys; 
  astro::SkyDir m_PulsarDir;
  CLHEP::Hep3Vector m_PulsarVectDir;
  
  const std::string& m_params; 

  //! Pulsar name  
  std::string m_PSRname;

  //! Name of file containing shape spectrum for model PSRShape
  std::string m_PSRShapeName;
  
  bool m_ff;

  //! Pulsar Coordinates and Directions
  double m_RA, m_dec;
  double m_l, m_b;  
  std::pair<double,double> m_GalDir;

  //! Ephemerides type
  std::string m_ephemType;

  //! Period and derivatives
  double m_period, m_pdot, m_p2dot;
  std::vector<double> m_periodVect, m_pdotVect, m_p2dotVect;

  //!Epoch and ephemerides validity range
  double m_t0, m_t0Init, m_t0End;
  std::vector<double> m_t0Vect, m_t0InitVect, m_t0EndVect;

  //! Frequency and derivatives
  double m_f0, m_f1, m_f2;
  std::vector<double> m_f0Vect, m_f1Vect, m_f2Vect;

  //! phase and turns at the epoch t0
  double m_phi0,m_N0;
  std::vector<double> m_phi0Vect,m_txbaryVect;

  //!Type of model
  int m_model;

  //!Pulsar flux
  double m_flux;

  //! model-dependent parameters
  int m_ppar0;
  double m_ppar1,m_ppar2,m_ppar3,m_ppar4;

  //!Minimum and maximum energy of the extracted photons
  double m_enphmin, m_enphmax;

  //! Random seed
  int m_seed;

  //! Flag for binary demodulation
  int m_BinaryFlag;

  //! Flag for enabling timing noise
  int m_TimingNoiseModel;

  //Binary parameters
  double m_Porb,m_asini,m_ecc,m_omega,m_t0PeriastrMJD,m_t0AscNodeMJD,m_PPN;

  //! Binary parameters - relative to PPN
  double m_Porb_dot,m_xdot,m_ecc_dot,m_omega_dot,m_gamma;
  double m_shapiro_r, m_shapiro_s;

};
#endif
