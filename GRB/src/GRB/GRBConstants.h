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
#include <iterator>
#include <iostream>
#include <fstream>

#include <math.h>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include "stdio.h"

#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TRandom.h"

#ifndef GRBCONSTANT_HH
#define GRBCONSTANT_HH 1

namespace cst
{
  static const double pi = acos(-1.0); 
  const Double_t c = 3.e+10; //cm
  const Double_t c2 = c*c;
  const Double_t mpc2  = 938.2;
  const Double_t mec2  = 0.510999;       //MeV
  const double st      = 6.65225e-25;
  const double gmin = 100;
  const double gmax = 10.*gmin;
  
  
  const float ab = .3; 
  const float ae = .3; 
  // const float ap = 1.-(ab-ae);
  
  const float csi = 1.0e-1;// Magnetic field fluctuation
  
  //  const double tmax = 10;
  
  const double emin =  1.0; //keV
  const double emax = 1e9;  //keV
  const double enph = 1e5;  //keV (100 MeV) 
  
  const    int Ebin =  100; 
  const    int Tbin =  100; 
  static const double de   = pow(emax/emin,1.0/Ebin);
  //  static const double dt   = tmax/(Tbin-1);
  const double erg2meV   = 624151.0;
  
  const double BATSE1=20.0;                   //20 keV
  const double BATSE2=1.0e+3;                 // 1 MeV
  const double GBM1=10.0;                     // 10 keV 
  const double GBM2=25.0e3;                   // 25 MeV 
  const double LAT1=50.0e3;                   // 50 MeV 
  const double LAT2=3.0e6;                    //300 GeV 
  //////////////////////////////////////////////////
};

class Parameters
{
 public:
  Parameters();
  double GetBATSEFluence();
  double GetLESI();
  double GetHESI();
  void SetNshell(int nshell);
  void SetFluence(double fluence);
  void SetEtot(double etot);
  void SetInitialSeparation(double initialSeparation);
  void SetInitialThickness(double initialThickness);		  
  //  void SetLESI(double beta);
  //void SetHESI(double alpha);
  void ReadParametersFromFile(std::string paramFile);
  void PrintParameters();
  inline long GetGRBNumber(){return m_GRBnumber;}
  inline void SetGRBNumber(long GRBnumber)
    {
      m_GRBnumber = GRBnumber;
      rnd->SetSeed(m_GRBnumber);
      rnd->Uniform();
    }
  
  TRandom *rnd;
  int    m_nshell;
  long m_GRBnumber;
  double m_fluence,m_etot,m_initialSeparation,m_initialThickness;
  //  double m_alpha,m_beta;
  
};

#endif
