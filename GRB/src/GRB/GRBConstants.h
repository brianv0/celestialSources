#ifndef GRBCONSTANT_HH
#define GRBCONSTANT_HH 1

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
  
  const double ab = 0.3;
  const double ae = 0.3; 

  const double emin =  10.0; //keV
  const double emax =  1e9;  //keV
  const double enph = 3.0e+4;  //keV (30 MeV) 
  
  const    int Ebin  =  50; 
  const double MinDT =  0.016; 
  static const double de   = pow(emax/emin,1.0/Ebin);
  const double erg2meV   = 624151.0;
  
  const double BATSE1=20.0;                   //20 keV
  const double BATSE2=50.0;;                 // 1 MeV
  const double BATSE3=100.0;                   //20 keV
  const double BATSE4=300.0;                 // 1 MeV
  const double BATSE5=1000.0;                 // 1 MeV
  const double GBM1=10.0;                     // 10 keV 
  const double GBM2=30.0e3;                   // 25 MeV 
  const double LAT1=30.0e3;                   // 50 MeV 
  const double LAT2=3.0e8;                    //300 GeV 
  //////////////////////////////////////////////////
};

class Parameters
{
 public:
  Parameters();
  ~Parameters(){ delete rnd;}
  int    GetNshell() {return m_Nshell;}
  double GetFluence(){return m_Fluence;}
  double GetEtot()   {return m_Etot;}
  double GetInitialSeparation(){return m_InitialSeparation;}
  double GetInitialThickness() {return m_InitialThickness;}
  double GetGammaMin(){return m_Gmin;}
  double GetGammaMax(){return m_Gmax;}
  double GetInverseCompton() {return m_InverseCompton;}

  double GetBATSEFluence();
  
  void SetGalDir(double l, double b);
  void SetNshell(int nshell);
  void SetFluence(double fluence);
  void SetEtot(double etot);
  void SetEpeak(double epeak);
  void SetInitialSeparation(double initialSeparation);
  void SetInitialThickness(double initialThickness);		  
  void SetGammaMin(double gmin);
  void SetGammaMax(double gmax);
  void SetInverseCompton(double ic);
  void ReadParametersFromFile(std::string paramFile, int NGRB=1);

  void ComputeParametersFromFile(std::string paramFile, int NGRB=1);

  void PrintParameters();
  inline UInt_t GetGRBNumber(){return m_GRBnumber;}
  inline std::pair<double,double> GetGalDir(){return m_GalDir;}
  inline double GetGamma(double gmin=0,double gmax=0)
    {
      if(gmin==0) gmin = m_Gmin;
      if(gmax==0) gmax = m_Gmax;
      return rnd->Uniform(gmin,gmax);
    }
  void SetGRBNumber(UInt_t GRBnumber);
  //  double GetNextPeak();
  //  inline   double GetTau(){return m_Tau;}

  TRandom *rnd; 

 private:

  UInt_t m_GRBnumber;
  int m_Type;
  int    m_Nshell;
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
