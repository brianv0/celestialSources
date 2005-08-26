#ifndef GRBCONSTANT_HH
#define GRBCONSTANT_HH 1

/*!   
 * \namespace cst 
 * \brief The namespace cst contains all the constant needed to the simulation.
 *
 * The parameter of the model are included in the GRBParam.txt and are computed by the Param class.
 *
 * \author Nicola Omodei       nicola.omodei@pi.infn.it 
 * \author Johann Cohen-Tanugi johann.cohen@pi.infn.it
 */
#include <vector> 
#include <string>
#include "TRandom.h"


namespace cst
{
  static const double pi = acos(-1.0); 
  static const Double_t c = 3.e+10; //cm
  static const Double_t c2 = c*c;
  static const Double_t mpc2  = 938.2;
  static const Double_t mec2  = 0.510999;       //MeV
  static const double st      = 6.65225e-25;
  static const double erg2meV   = 624151.0;

  /// Fraction of energy that goes in magnetic field  
  const double ab  = 0.3;
  /// Fraction of energy that goes in accelerating electrons  
  const double ae  = 0.3; 
  /// Fraction of electron which are accelerated by the shock
  const double csi = 0.1; 
  /// minimum energy for the computed spectrum
  const double emin =  10.0; //keV
  /// maximum energy
  const double emax =  1e9;  //keV
  /// default value for the minimum energy of the extracted photons (LAT photons)
  const double enph = 3.0e+4;  //keV (30 MeV) 
  /// Number of energy bin (logarithically spaced)
  const    int Ebin  =  50; 
  /// Time resolution for the GBM spectra.
/// Temporal resolution
  const    double  TimeBinWidth   =  0.016; //s 1 msec
  /// Time resolution for the GBM spectra.
  const    double  GBMTimeBinWidth   =  0.016; //s 16 msec
  static const double de   = pow(emax/emin,1.0/Ebin);
  /// Bottom edge of the 1st channel of CGRO/BATSE 20 keV
  const double BATSE1=20.0;       
  /// Top edge of the 1st channel of CGRO/BATSE 50 keV
  const double BATSE2=50.0;                 
  /// Top edge of the 2nd channel of CGRO/BATSE 100 keV
  const double BATSE3=100.0;                
  /// Top edge of the 3rd channel of CGRO/BATSE 300 keV
  const double BATSE4=300.0;                
  /// Top edge of the 4th channel of CGRO/BATSE 1MeV 
  const double BATSE5=1000.0;                
  /// Low edge of GLAST/GBM energy band 10 keV 
  const double GBM1=10.0;                     
  /// Top edge of the of GLAST/GBM energy band 30 MeV
  const double GBM2=30.0e3;                   
  /// Low edge of GLAST/LAT energy band 30 MeV 
  const double LAT1=30.0e3;                   
  /// Top edge of the of GLAST/LAT energy band 300 GeV
  const double LAT2=3.0e8;                    
};
/*!
 * \class Parameters
 * \brief This class reads the parameters from the file, compute the parameters of the physical model 
 * 
 */
class Parameters
{
 public:
  /// Constructor: initializes the (root) random number generator
  Parameters();
  /// Destructor: delete the random number generator
  ~Parameters(){ delete rnd;}
  /*!
    The parameters file Param.txt contains the following observable quantities (typical values between parenthesis):
    - Fluence between 20 keV and 1 MeV, in erg/cm^2 (BATSE values)
    - Burst type: S=short bursts, L = long burst, R -> 25% S, 75% L
    - Cut-off energy \f$ E_{co} \f$, in GeV. (3, 10 GeV) 
    - Peak Energy \f$ E_{peak} \f$ of the synchrotron spectrum, in keV 
    (it is sampled from a log-normal distribution with mean 235 keV and sigma 1.75. Moreover, \f$ E_{peak,SHORT}=E_{peak,LONG}/2 \f$)
    - IC/Syn is the ratio between the Inverse Compton peak and the Synchrotron peak of the 
    \f$ e^2 N(e)\f$ spectrum (0=pure synchrotron. < 10 typically)
    - GBM flag, if 1 the program generates a series of files to be used with the GBM software
    
    The variability scale \f$ t_v \f$ is set by the observations, requiring that for short bursts it is approximately equal to the burst duration 
    (\f$ t_v\approx T_{90}\f$), whicle for long bursts it is sampled from the observed distribution of the FWHM: (\f$ t_v=10^{Gaus(0.0,0.5)}\f$).
    The variability time scale (which is directly related to the pulse width) is contsrained to a minimum: 
    The shortest variability time scale is a fraction (1/50) of the duration of the bursts: 
    long bursts will likely have long pulses instead of many short pulses.
    
    This quantities are converted into the parameters of the model by a series of relations:

    - \f$ \Gamma=40.5 * E_{co} \f$
    - \f$ \frac{\Delta_0}{10^7}  =  13.6~(3~\alpha_B)(\frac{E_{tot}}{10^{52}})(\frac{3 \alpha_e \xi}{E_{co}})(\frac{E_{peak}}{100} t_v)^2 \f$
    - \f$ r_0  =  2~c~t_v \f$
    
    The total energy is set (and fixed) to  \f$ 10^{52} \f$ erg, and \f$ \Gamma_M/\Gamma_m=2 \f$, so that: 
    \f$ \Gamma_m = \frac{2~\Gamma}{\Gamma_M/\Gamma_m+1} \f$, and \f$ \Gamma_M= \frac{\Gamma_M}{\Gamma_m}\Gamma_m \f$ 
    This is a model dependent way to correlate observables and model parameters.
  */
  void ComputeParametersFromFile(std::string paramFile, int NGRB=1);
  /// Print out method
  void PrintParameters();
  /// Return a boolean true is the generation of the GBM output has been chosen
  inline bool GenerateGBM(){return m_GBM;}
  /// Return a GRB number identifier, used to initialize the random number generator
  inline UInt_t GetGRBNumber(){return m_GRBnumber;}
  /// return galacic l and b coordinates of the burst
  inline std::pair<double,double> GetGalDir(){return m_GalDir;}
  /// return The fluence in erg/cm^2 between 20 keV and 1 MeV
  inline double GetFluence(){return m_Fluence;}
  /// return the duration
  inline double GetDuration(){return m_Duration;}
  /// return the intrinsic energy (energy of the source) expressed erg
  inline double GetEtot()   {return m_Etot;}
  /// return the initial separation between shells as they are emitted from the central engine. (cm)
  inline double GetInitialSeparation(){return m_InitialSeparation;}
  /// Initial thickness of the shells (cm)
  inline double GetInitialThickness() {return m_InitialThickness;}
  /// Minimum lorentz factor 
  inline double GetGammaMin(){return m_Gmin;}
  /// Maximum lorentz factor 
  inline double GetGammaMax(){return m_Gmax;}
  /*!
    Returns the Inverse Compton parameter.
  */
  inline double GetInverseCompton() {return m_InverseCompton;}
  /*! 
    Return a value of the fluence sampled from the BATSE distribution. 
    Different distributions are considered for long and short bursts.
    
    \f$ 10^{Gaus(-6.3,0.57)} \f$ erg/cm^2 (Short Bursts)
    \f$ 10^{Gaus(-5.4,0.62)} \f$ erg/cm^2 (Long Bursts)
  */
  double GetBATSEFluence();
  /*! 
    Return a value of the duration sampled from the BATSE T90 distribution. 
    Different distributions are considered for long and short bursts.
    
    \f$ 10^{Gaus(-0.2,0.55)} \f$ s (Short Bursts)
    \f$ 10^{Gaus(1.46,0.49)} \f$ s (Long Bursts)
    
  */
  double GetBATSEDuration();
  
  void SetGalDir(double l, double b);
  void SetFluence(double fluence);
  void SetEtot(double etot);
  void SetInitialSeparation(double initialSeparation);
  void SetInitialThickness(double initialThickness);		  
  void SetGammaMin(double gmin);
  void SetGammaMax(double gmax);
  void SetInverseCompton(double ic);
  inline void SetGBMOutput(bool flag){m_GBM=flag;}
  inline void SetQGOutput(bool flag){m_QG=flag;}
  /// Initialize the random number generator
  void SetGRBNumber(UInt_t GRBnumber);
  
  /*!
    This method sets the basic parameters of the model:
    They are:
    -# The fluence of the burst
    -# The total energy (\f$ E_{tot} \f$)
    -# The initial radius of the shells (\f$ r_0 \f$)
    -# The initial thickness (\f$ \Delta_0 \f$) 
    -# The minimum Lorentz Factor of the shells  (\f$\Gamma_m \f$)
    -# The maximum Lorentz Factor of the shells  (\f$\Gamma_M \f$)
    -# The Inverse Compton parameter (ratio between the IC peak and the synchrotron peak of the \f$e^2N(e)\f$ flux)
  */
  void SetParameters(double fluence, 
		     double etot, 
		     double r0, 
		     double dr0, 
		     double gmin, 
		     double gmax, 
		     double ic);
  /// The random number generator.
  TRandom *rnd; 

  const bool QG() const {return m_QG;}


 private:
  
  bool m_GBM;
  bool m_QG;
  UInt_t m_GRBnumber;
  int m_Type;
  double m_Duration;
  double m_Gmin;
  double m_Gmax ;
  double m_Fluence;
  double m_Ep;
  double m_Etot   ;
  double m_InitialSeparation;
  double m_InitialThickness ;
  double m_InverseCompton ;
  //  double m_Tau ;
  std::pair<double,double> m_GalDir;
};

#endif
