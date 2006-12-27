#ifndef GRBOBSCONSTANT_HH
#define GRBOBSCONSTANT_HH 1

#include <vector> 
#include <string>
#include "TRandom.h"
#include "TF1.h"

/*! 
 * \namespace ObsCst
 * \brief Namespace containing the constants of the model such as the binning in time and energy
 
 * \author Nicola Omodei       nicola.omodei@pi.infn.it 
 */
namespace ObsCst
{
  /// minimum energy for the computed spectrum (keV)
  const double emin = 10.0;
  /// maximum energy (keV)
  const double emax = 3e8;
  /// default value for the minimum energy of the extracted photons -LAT photons- (keV)
  const double enph = 3e4; 
  //////////////////////////////////////////////////
  ///reference energy for the Universal Width (keV)
  const double E0           = 20.0;  
  /// power law  index for the full width at half maximum - energy relation: \f$ W(e)=W0 (e/E0)^{-We} \f$
  const double We           = 0.33;//0.33; 
  /// for the displacements between peaks
  const double deltaTPeak   = 0.0;//0.5;  //for the displacements between peaks
  /// Number of energy bins (logarithmically spaced)
  const    int Ebin =  50; 
  /// Time resolution for the GBM spectra.
  const    double  TimeBinWidth   =  0.016; //s
  static const double de   = pow(emax/emin,1.0/Ebin);
  
  const double erg2meV   = 624151.0;
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
  //////////////////////////////////////////////////
};

/*! 
  \class GRBobsParameters
  Class instantiated to access general parameters and constants.
  
  This classes uses observed quantities to generate the parameters for the phenomenological model.
 
 
 * \author Nicola Omodei       nicola.omodei@pi.infn.it 
 */


class GRBobsParameters
{
 public:
  /// Constructor. Initialized the random number generator.
  GRBobsParameters();
  /// Destructor
  ~GRBobsParameters()
    { 
      delete rnd;
    }

  inline  double GetFluence()        {return m_fluence;}
  inline  double GetRiseTime()       {return m_riseTime;}
  inline  double GetDecayTime()      {return m_decayTime;}
  inline  double GetPulseHeight()    {return m_pulseHeight;}
  inline  double GetPulseSeparation(){return m_pulseSeparation;}
  inline  double GetPeakedness()     {return m_Peakedness;}
  inline  double GetEpeak()          {return m_Epeak;}
  inline  double GetLowEnergy()      {return m_LowEnergy;}
  inline  double GetHighEnergy()     {return m_HighEnergy;}
  inline  long   GetGRBNumber()      {return m_GRBnumber;}
  inline  double GetDuration()       {return m_duration;}
  
  inline std::pair<double,double> GetGalDir(){return m_GalDir;}
  //////////////////////////////////////////////////
  void   SetGRBNumber(long);
  void   SetFluence(double);
  void   SetDuration(double);
  /// Sets the low energy spectral index and the high energy spectral index.
  /// \param alpha low energy spectral index (\f$-3<\alpha<1\f$)
  /// \param beta high energy spectral index (\f$\beta>-1\f$,\f$\beta>\alpha\f$)
  /// If alpha is not in the interval than its value is sampled from a gaussian distribution centered in -1 with sigma 0.4 
  /// If beta is not in the interval than its value is sampled from a gaussian distribution centered in -2.25 with sigma 0.4 
  /// Notice that alpha and beta are scaled with the factor We, so that they represent the real
  /// measured values, once the spectral variation introduced by the energy dependent temporal profile is taken into account.
  inline void SetAlphaBeta(double alpha, double beta)
    {
      m_LowEnergy       = alpha;
      m_HighEnergy      = beta;
      
      while (m_LowEnergy<-3.0 || m_LowEnergy>1.0)                        m_LowEnergy  = rnd->Gaus(-1.0,0.4);
      while(m_HighEnergy >= m_LowEnergy || m_HighEnergy >= -1.0)         m_HighEnergy = rnd->Gaus(-2.25,0.4);
      
      m_LowEnergy       += ObsCst::We;
      m_HighEnergy      += ObsCst::We;
    } 
  /// Set the minimum photon energy for generating LAT photons.
  void   SetMinPhotonEnergy(double);
  /// Set the galactic position in the sky of the GRB in <em>(l,b)</em> coordinates.
  void   SetGalDir(double, double);  
  /// The distribution of fluences is taken fron the BATSE catalog.
  /// \f$ N(\log(f_{tot}))=A \exp[-0.5\left(\frac{\log(f_{tot})-\overline{\log(f_{tot})}}{\sigma_{f_{tot}}}\right)^2]\f$ <br>
  /// Where, for short bursts: \f$\overline{\log(f_{tot})}= -0.20\f$, \f$\sigma_{f_{tot}}=0.55\f$.
  ///        for long  bursts: \f$\overline{\log(f_{tot})}=  1.46\f$, \f$\sigma_{f_{tot}}=0.49\f$.
  ///
  double GetBATSEFluence();

  /// Distribution of the bursts duration expressed by the quantity \f$T_{90}\f$, taken from the BATSE catalog. 
  /// \f$N(\log(T_{90}))=A \exp[-0.5\left(\frac{\log(T_{90})-\overline{\log( T_{90})}}{\sigma_{T_{90}}}\right)^2]\f$ <br>
  /// Where, for short bursts: \f$\overline{\log(T_{90})}= -6.35\f$, \f$\sigma_{T_{90}}=0.57\f$.
  ///        for long  bursts: \f$\overline{\log(T_{90})}= -5.39\f$, \f$\sigma_{T_{90}}=0.62\f$.
  double GetBATSEDuration();  
  
  /*!
    When it is called, It random generates a set of parameter, to be used with a GRBobsPulse.
    - The peakedness is set to 1.5
    - The FWHM at 20 keV is: for short burst is a fraction of its duration (\f$T_{90}\f$): 
    \f$\xi\times T_{90}\f$ where \f$\xi\f$ is a (uniform) random number between 0 and 1.
    For long bursts the pulse duration is sampled from a log normal distruibution: <em> pow(10.0,rnd->Gaus(-0.1,0.5))</em>
    - The rise to decaty ratio (\f$\rho\f$) is randomly sampled from a gaussian of mean 0.4 and sigma 0.1 (<em> rnd->Gaus(0.4,0.1))</em>).
    - The rise time is then: \f$\rho/(1+\rho)*(log(2))^{-1/\nu}*FWHM\f$;
    - The decay time is: \f$1/(1+\rho)*(log(2))^{-1/\nu}*FWHM\f$;
    - The Pulse Height is uniformly sampled (<em> rnd->Uniform() </em>).
    - The Peak energy is sampled from a log-normal distribution <em> pow(10.,rnd->Gaus(log10(235.0),log10(1.75)))</em>. 
    If the burst is long the peak energy is shifted to lower energies.
  */
  void GenerateParameters();
  
  /// Read parameters from the parameter file. Each line correspond to a GRB parameters set. 
  /// This method is called only in the test program GRBROOTTest.cxx, while GRBobsmanager read the parameter parsing the xml string.
  void ReadParametersFromFile(std::string paramFile, int NGRB=1);

  //// Print out utility for debugging.
  void PrintParameters();

  TRandom *rnd;
 
 private:
  int m_Type;
  double  m_Peakedness;
  double  m_FWHM;
  double m_pulseSeparation;

  double  m_Epeak;
  double m_LowEnergy,m_HighEnergy;
  double m_fluence;
  double m_duration;
  double m_RD,m_riseTime,m_decayTime;
  double m_pulseHeight;

  long   m_GRBnumber;
  double m_enph;
  std::pair<double,double> m_GalDir;
};

#endif
