/*! 
 * \class GRBConstants
 * \brief Class instantiated to access general parameters and constants.
 *  
 * The namespace cst contains all the constant needed to the simulation.
 * All the constant relative to physical model are included in the GRBParam.txt
 *
 * \author Nicola Omodei       nicola.omodei@pi.infn.it 
 * \author Johann Cohen-Tanugi johann.cohen@pi.infn.it
 */

#ifndef GRBCONSTANTS_HH
#define GRBCONSTANTS_HH 1
#include <string>
#include <cmath>
#include "CLHEP/Random/RandomEngine.h"

/*!
  \brief Namespace contaning all the constants
  
  In this namespace are stored all the constants needed by the GRB Physical 
  simulator.
  The constants contained here usually should not be changed by a general user.
  
*/
namespace cst
{
  // Universal constants :
  //! Proton Rest Mass (MeV)
  const double mpc2      = 938.2; 
  //! Electron Rest Mass (MeV)
  const double mec2      = 0.510999;       //MeV
  //! Conversion fron \e erg to \e MeV
  const double erg2MeV   = 624151.0;
  //! Conversion fron \e Mpc to \e cm
  const double mpc2cm    = 3.0857e+24;
  // Conversion fron \e Gauss to \e MeV
  //  const double G2MeV     = 815.78;
  
  //! Thompson cross section (\f$ cm^2 \f$)
  const double st        = 6.65225e-25;
  //! Magnetic permeability in \f$ N A^{-2}\f$
  const double mu0       = 4.0*M_PI*1.0e-7; 
  //! Reference value for the magnetic field \e Gauss. 
  const double   BQ        = 4.413e+13; 
  //! Speed of light in \f$ cm/sec\f$
  const double c         = 2.98e+10;
  //! Square of the speed of light
  const double c2        = c*c;
  //! Planck constant in eV*sec.
  const double hplanck   = 4.13567e-15;
  //! \f$\pi\f$
  const double pi        = 3.1415926535897932385;
  //! Hobble`s constant in \f$km/(s Mpc)\f$
  const double Hubble    = 6.5e+1; 
  //! w Zeldovich
  const double wzel      = 0.0;
  //! Scale factor, indicates the ratio of accelerated electrons
  const double csi       = 1.0; 
  //! Part of the internal energies that goes in electrons
  const double alphae    = .5;
  //! Part of the internal energies that goes in magnetic fiels
  const double alphab    = .1;
  //! The shock accelerate electron with a power low energy disptribution:
  //! \f$N(E)dE\propto E^{-p}dE\f$ 
  const double p         = 2.5;
  //! Viscosity of the circum burst material
  const double viscosity = 0.0;
  //! Number of time steps
  const int    nstep        = 300; 
  //! Number of energy steps
  const int    enstep       =  50;
  //! Minimum energy at which the simulator will compute the flux 
  const double enmin     = 1.0e+3;
  //! Maximum energy at which the simulator will compute the flux 
  const double enmax     =1.0e+12;
  // Minimal Temporal separation between 2 photons
  // const double DeadTime  = 1.0e-5; //sec
  
  /*! Compton Scattering calculation

    This constant give an estimation of how Compton scattering will 
    partecipate to the spectrum. flagIC has to be in the interval [0,1],
    if ==0, no inverse compton will be calculated. If == 1 the IC effect is
    calculated tacking into account the optical depht of the shell.
    \sa GRBSynchrotron, GRBICompton, GRBSim
  */
  const float flagIC     = 1;
  /*! \brief Flag to compute the Quantum Gravity Effect
     
     Quantum gravity, if present, will calculate the dispersion low 
     on the arrival time for the photons depending on thei energies.
     \sa RadiationProcess::timeShiftForDispersion()
     
   */
  const bool  flagQG     = false;
  /*!Flag for saving in output file
    
    If true it save into an output file. It is needed by 
    <a href="../../src/test/other/GRBROOTtest.cxx"> GRBROOTtest.cxx</a>
  */
  const bool savef=false;
  //! Name of the output file
  const std::string paramFile= "GRBdata.txt";
  /*! Indicates the pulse shape

    Different pulse shape are considered.
    The pulse shape depends on how the particle are accelerated in a shock 
    as function of time.
    - sgauss = Simmetric pulse shape. It is a double exponential with 
    the two characteristic time equal to the cooling time of the radiative
     process. No intrinsic delay (or 'lag' has been considered between 
    different energies.
    - agauss = asimmetric exponential function, the rise time and the decay
    time are different. The peak time depend on the energy, 
    and this give up a delay between high energies and low energies. 
    -else . any other choiche set up a default function expressed by a FRED
     like function. A Fast Rise pulse followed by an Exponential Decay.
  \sa RadiationProcess::electronNumber
  */
  const std::string pulse_shape="agauss"; 
}

class GRBConstants 
{ 
 public:

  GRBConstants();

  ~GRBConstants() { }

  //! Initialize the Random Number Generator. Every run is differrent.
  static void InitializeRandom(long seed = 0);
  
  static HepRandomEngine* GetTheRandomEngine(long seed = 0);
  
  //! Parameters are read from a file using facilities::Util::expandEnvVar method
  int ReadParam();
  
  //! Printout Utility
  void Print();
  
  //! Save the parameters in a file. It could be usefull to mantain trace of the models runned.
  void Save(bool flag=false);


  /*!  Defines the engine`s type

    The type of engine defines the way GRBengine creates the shocks vector

    \param value define the type: 
    - type 0: Observed parameters
    - type 1: Reads the Physical parameters and create a shock with 
    the defined parameters
    - type 2: Reads the physical parameters and generates shocks 
    starting from the collision of two shells of the specified parameters.

    \sa GRBengine, GRBShock
   */
  inline void setEngineType(int value){m_engineType = value;}
  
  //! Retrn The engine type 
  inline int EngineType(){return m_engineType;}
    
  ////////////////////  Engine: Shell Generator

  //  inline int Nshell() {return m_nshell;}

  //  inline void setNshell(int value=10){m_nshell = value;}
    
  ////////////////////  Engine: Shock Generator

  inline int Nshock() {return m_nshock;}

  inline void setNshock(int value=10){m_nshock = value;}

  inline double Duration(){return m_duration;}

  inline void setDuration(double value=0.4){m_duration = value;}

  inline double RiseTime(){return m_rt;}

  inline void setRiseTime(double value=0.4){m_rt = value;}
 
  inline double DecayTime(){return m_dt;}
 
  inline void setDecayTime(double value=0.4){m_dt = value;}
 
  inline double PeakEnergy(){return m_peak;}
 
  inline void setPeakEnergy(double value=0.4){m_peak = value;}
  /*! Sets the type of shell
    
    GRBShell can have different geometries. 
    This function set the appropriate geometry
    - Type 0: Isotropic emission
    - Type 1: Jet emission
  */
  inline void setShellType(int value){m_shellType = value;}
  
  //! Returns the shell`s type
  inline int ShellType(){return m_shellType;}
  
  //////////////////// Shell: Spherical Shells

  inline double ShellRadius(){return m_d0;}
  
  inline void setShellRadius(double value){m_d0 = value;}

  //////////////////// Shell: Jet Shells
 
  inline double JetRadius(){return m_r0;}
  
  inline void setJetRadius(double value){m_r0 = value;}
  
  inline double JetAngle(){return m_angle;}
  
  inline void setJetAngle(double value){m_angle = M_PI/180. * value;}
  
  //////////////////// Common
 
  inline double Thickness(){return m_t0;}

  inline void setThickness(double value){m_t0 = value;}
   
  inline double Etot(){return m_etot;}

  inline void setEtot(double value){m_etot = value;}

  inline double GammaMin(){return m_g0;}

  inline void setGammaMin(double value=100.0){m_g0 = value;}

  inline double GammaMax(){return m_g1;}

  inline void setGammaMax(double value=1000.0){m_g1 = value;}

  inline double Redshift() {return m_redshift;}

  inline void setRedshift(double value=1.0){m_redshift = value;}
  /*! Set the minimum energy of the extracted photons 
    
    All the photons drawn from the spectrum will be energy greater then m_enph
  */
  inline void setEnergyPh(double value=25.0e+3){m_enph = value;}
  
  inline double EnergyPh(){return m_enph;}
  
  double Distance();
  
  inline std::string GetParamFile() {return  m_paramFile;}
 
 private:
  char* m_burst_type;
  //  int m_nshell;
  double m_d0;
  int m_nshock;
  double m_r0;
  double m_angle;
  double m_redshift;
  double m_etot;
  double m_g0,m_g1;
  double m_t0;
  double m_duration;
  double m_rt,m_dt,m_peak;
  double m_enph;
  std::string m_paramFile;
  int m_engineType;
  int m_shellType;
};

#endif


