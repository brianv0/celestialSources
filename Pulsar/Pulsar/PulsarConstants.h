#ifndef PulsarCONSTANT_HH
#define PulsarCONSTANT_HH 1

/*! 
 * \class PulsarConstants
 * \brief Class instantiated to access general parameters and constants.
 *  
 * The namespace cst contains all the constant needed to the simulation.
 *
 * \author Nicola Omodei       nicola.omodei@pi.infn.it 
 * \author Massimiliano Razzano massimilia.razzano@pi.infn.it
 */

#include <vector> 
#include <string>
#include "TRandom.h"


namespace cst
{
  const double EnNormMin = 1e5; // KeV (100 MeV) //Normalization Interval
  const double EnNormMax = 1e7; // KeV (10 GeV)
  
  const    int Ebin =  50 ; // Energy bins 
  const    int Tbin =  200; // Time bins 

  const double erg2meV   = 624151.0;


  const double GBM1=10.0;                     // 10 keV 
  const double GBM2=25.0e3;                   // 25 MeV 
  const double LAT1=30.0e3;                   // 30 MeV 
  const double LAT2=300.0e6;                  // 300 GeV 
  const double EGRET1=30.0e3;                 // 30 MeV
  const double EGRET2=100.0e3;                // 100MeV
  const double EGRET3=10.0e6;                 // 10 GeV

  const double StartMissionDateMJD = 51910.0; //Start Date: 1 Jan 2001, 00:00.00
  const double JDminusMJD = 2400000.5; //Difference between JD and MJD
  const int SecsOneDay = 86400; //seconds in one day

  //////////////////////////////////////////////////
  
};

#endif
