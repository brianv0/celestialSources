/*! 
 * \class GRBConstants
 * \brief Class instantiated to access general parameters and constants.
 *  
 * The namespace cst contains all the constant needed to the simulation.
 * All the constant relative to physical model are included in the GRBParam.txt
 *
 * \author Nicola Omodei nicola.omodei@pi.infn.it
 * \author Johann Cohen-Tanugi johann.cohen@pi.infn.it
 */

#ifndef GRBCONSTANTS_HH
#define GRBCONSTANTS_HH 1
#include <string>
#include <cmath>
#include "CLHEP/Random/RandomEngine.h"
namespace cst
{
  // Universal constants :
  // Electron Rest Mass (MeV)
  const double mpc2      = 938.2; 
  const double erg2MeV   = 624151.0;
  const double mpc2cm    = 3.0857e+24;

  const double G2MeV     = 815.78;
  const double st        = 6.65225e-25;
  const double mec2      = 0.510999;       //MeV
  const double mu0       = 4.0*M_PI*1.0e-7; //N A^-2
  const double   BQ        = 4.413e+13; // Gauss
  // light speed in cm/sec.
  const double c         = 2.98e+10;
  const double c2        = c*c;
  // Planck constant in eV*sec.
  const double hplanck   = 4.13567e-15;
  const double pi        = 3.1415926535897932385;
  const double Hubble    = 6.5e+1; 
  const double wzel      = 0.0;

  const double csi       = 1.0; //csi = 1, more efficient IC
  const double alphae    = .5;
  const double alphab    = .1; //smaller is alphab greater is the IC efficiency
  const double p         = 2.5;
  const double viscosity = 0.0;
  // Internal Parameters
  const int    nstep        = 400;
  const int    enstep       = 100;
  const double enmin     = 1.0e+3;
  const double enmax     =1.0e+12;
  // Minimal Temporal separation between 2 photons
  const double DeadTime  = 1.0e-5; //sec
  // flag =[0,1], if ==0, No inverse compton;
  const float flagIC     = 1;
  // Flag to compute the Quantum Gravity Effect
  const bool  flagQG     = false;
  // If true it save into an output file 
  const bool savef=false;
  // Name of the output file
  const std::string paramFile= "GRBdata.txt";
  /* Different pulse shape are considered.
   * The pulse shape depends on how the particle are accelerated in a shock 
   * as function of time.
   * - sgauss = Simmetric pulse shape. It is a double exponential with 
   * the two characteristic time equal to the cooling time of the radiative
   *  process. No intrinsic delay (or 'lag' has been considered between 
   * different energies.
   * - agauss = asimmetric exponential function, the rise time and the decay
   * time are different. The peak time depend on the energy, 
   * and this give up a delay between high energies and low energies. 
   * -else . any other choiche set up a default function expressed by a FRED
   *  like function. A Fast Rise pulse followed by an Exponential Decay.
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

  //! Defines the engine`s type
  inline int EngineType(){return m_engineType;}
  
  inline void setEngineType(int value){m_engineType = value;}
  
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
  //! Defines the shell`s type
  inline int ShellType(){return m_shellType;}
  
  inline void setShellType(int value){m_shellType = value;}
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
  //! Set the minimum energy of the extracted photons 
  //! (all the photons drown from the spectrum will be energy greater then m_enph)
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


