#ifndef GRBOBSCONSTANT_HH
#define GRBOBSCONSTANT_HH 1

#include <vector> 
#include <string>
#include "TRandom.h"
#include "TF1.h"

/*! 
 * \namespace ObsCst
 * \brief Namespace containing the constants of the model such as the binning in tiime and energy
 
 * \author Nicola Omodei       nicola.omodei@pi.infn.it 
 */
namespace ObsCst
{
  const double emin = 10.0; //keV
  const double emax = 1e8;  //keV
  const double enph = 1e5;  //keV (100 MeV) 
  
  const    int Ebin =  50; 
  const    double  TimeBinWidth   =  0.016; //s
  static const double de   = pow(emax/emin,1.0/Ebin);

  const double erg2meV   = 624151.0;
  
  const double BATSE1=20.0;                   //20 keV
  const double BATSE2=50.0;;                 // 1 MeV
  const double BATSE3=100.0;                   //20 keV
  const double BATSE4=300.0;                 // 1 MeV
  const double BATSE5=1000.0;                 // 1 MeV
  const double GBM1=10.0;                     // 10 keV 
  const double GBM2=30.0e3;                   // 25 MeV 
  const double LAT1=50.0e3;                   // 50 MeV 
  const double LAT2=3.0e6;                    //300 GeV 
  //////////////////////////////////////////////////
};

/*! 
 * \class GRBobsParameters
 * \brief Class instantiated to access general parameters and constants.
 
 * \author Nicola Omodei       nicola.omodei@pi.infn.it 
 */


class GRBobsParameters
{
 public:
  
  GRBobsParameters();
  ~GRBobsParameters()
    { 
      delete rnd;
      delete PeakSeparationDist;
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
  inline  int    GetNumberOfPulses() {return m_numberOfPulses;}
  
  inline std::pair<double,double> GetGalDir(){return m_GalDir;}
  //////////////////////////////////////////////////
  void   SetGRBNumber(long);
  void   SetFluence(double);
  void   SetRiseTime(double);
  void   SetDecayTime(double);
  void   SetPulseHeight(double);
  void   SetPulseSeparation(double);
  void   SetNumberOfPulses(int);
  void   SetTau(double);
  void   SetMinPhotonEnergy(double);
  void   SetGalDir(double, double);  
  double GetBATSEFluence();

  void GenerateParameters();

  void ReadParametersFromFile(std::string paramFile, int NGRB=1);

  void PrintParameters();

  TRandom *rnd;
 
 private:
  double  m_Peakedness;
  double  m_FWHM;
  double  m_Epeak;
  double m_LowEnergy,m_HighEnergy;
  double m_fluence;
  int m_numberOfPulses;
  double m_riseTime;
  double m_decayTime;
  double m_pulseHeight;
  double m_pulseSeparation;
  long   m_GRBnumber;
  double m_Tau;
  double m_enph;
  std::pair<double,double> m_GalDir;
  TF1 *PeakSeparationDist;
};

#endif
