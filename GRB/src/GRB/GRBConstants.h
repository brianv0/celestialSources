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
  const double alphae    = .33;
  const double alphab    = .33; //smaller is alphab greater is the IC efficiency
  const double p         = 2.5;
  const double viscosity = 0.0;
  // Internal Parameters
  const int    nstep        = 400;
  const int    enstep       = 50;
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
}

class GRBConstants 
{ 
 public:
  //! Constructor calls method readParam wich gets an external file
  GRBConstants();

  //! Destructor
  ~GRBConstants() { }

  //! Initialize the Random Number Generator. Every run is differrent.
  static void InitializeRandom(long seed = 0);
  
  //! Returns the random engine
  static HepRandomEngine* GetTheRandomEngine(long seed = 0);
  
  //! Parameters are read from a file using facilities::Util::expandEnvVar method
  int ReadParam();
  
  //! Printout Utility
  void Print();
  
  //! Save the parameters in a file. It could be usefull to mantain trace of the models runned.
  void Save(bool flag=false);
  
  ////////////////////  Engine: Shell Generator
  //! Number of shells generated from the source
  inline int Nshell() {return m_nshell;}
  //! Set the number of shells.
  inline void setNshell(int value=10){m_nshell = value;}
  //! Initial separation between shells (cm)
  inline double ShellRadius(){return m_d0;}
  //! Set the initial separation between the shells.
  inline void setShellRadius(double value){m_d0 = value;}
    
  ////////////////////  Engine: Shock Generator
  //! Number of shells generated from the source
  inline int Nshock() {return m_nshock;}
  //! Set the number of shells. If the arguments is 0 it returns a random number.
  inline void setNshock(int value=10){m_nshock = value;}
  //! Duration of the burst
  inline double Duration(){return m_duration;}
  //! Set the duration of the burst
  inline void setDuration(double value=0.4){m_duration = value;}
  //! Duration of the burst
  inline double RiseTime(){return m_rt;}
  //! Set the duration of the burst
  inline void setRiseTime(double value=0.4){m_rt = value;}
  //! Duration of the burst
  inline double DecayTime(){return m_dt;}
  //! Set the duration of the burst
  inline void setDecayTime(double value=0.4){m_dt = value;}
  //! Duration of the burst
  inline double PeakEnergy(){return m_peak;}
  //! Set the duration of the burst
  inline void setPeakEnergy(double value=0.4){m_peak = value;}
  
  //////////////////// Shell: Spherical Shells
  //////////////////// Shell: Jet Shells
  //! Radius of the Jet (cm)
  inline double JetRadius(){return m_r0;}
  //! Sets the Radius of the Jet (cm)
  inline void setJetRadius(double value){m_r0 = value;}
  //! Returns the angle between the jet direction and the angle of sight, in radiants.
  inline double JetAngle(){return m_angle;}
  //! Sets the jet angle, \param value is in degree.
  inline void setJetAngle(double value){m_angle = M_PI/180. * value;}
  
  //////////////////// Common
  //! Thickness of the shells
  inline double Thickness(){return m_t0;}
  //! Set the initial thickness
  inline void setThickness(double value){m_t0 = value;}
   
  //! Total Energy available (ergs)
  inline double Etot(){return m_etot;}
  //! Set the total Energy available.
  inline void setEtot(double value){m_etot = value;}
  //! Minimum Lorentz factor of the shells.
  inline double GammaMin(){return m_g0;}
  //! Set the minimum Lorentz factor of the shells
  inline void setGammaMin(double value=100.0){m_g0 = value;}
  //! Maximum Lorentz of the shells
  inline double GammaMax(){return m_g1;}
  //! Set the maximum Lorenz of the shells
  inline void setGammaMax(double value=1000.0){m_g1 = value;}

  //! Redshift of the source.
  inline double Redshift() {return m_redshift;}
  //! Set the redshift of the source. If the arguments is 0 it returns a random number.
  inline void setRedshift(double value=1.0){m_redshift = value;}
  //! Set the minimum energy of the extracted photons (all the photons drown from the spectrum will be energy greater then m_enph)
  inline void setEnergyPh(double value=25.0e+3){m_enph = value;}
  //! Return the minimum energy of the extracted photons.
  inline double EnergyPh(){return m_enph;}
  //! Return the distance of a GRB from the Earth
  double Distance();
  //! Return the name of the file where the partameter will be saved.
  inline std::string GetParamFile() {return  m_paramFile;}
 
 private:
  char* m_burst_type;
  int m_nshell;
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
};

#endif


