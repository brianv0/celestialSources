#ifndef PulsarCONSTANT_HH
#define PulsarCONSTANT_HH 1

/*! 
 * \class PulsarConstants
 * \brief Class instantiated to access general parameters and constants.
 *  
 * The namespace cst contains all the constant needed to the simulation.
 * All the constant relative to physical model are included in the PulsarParam.txt
 *
 * \author Nicola Omodei       nicola.omodei@pi.infn.it 
 * \author Johann Cohen-Tanugi johann.cohen@pi.infn.it
 */
//#include <iterator>
//#include <iostream.h>
//#include <fstream.h>

/* #include <math.h> */
#include <vector> 
#include <string>
/* #include <algorithm> */
/* #include <cmath> */
/* #include "stdio.h" */

/* #include "TFile.h" */
/* #include "TCanvas.h" */
/* #include "TPad.h" */
/* #include "TF1.h" */
/* #include "TF2.h" */
/* #include "TH1D.h" */
/* #include "TH2D.h" e*/
#include "TRandom.h"


namespace cst
{
  const double emin = 1e5; //KeV (100 MeV)
  const double emax = 1e8;  //KeV (100 GeV)
  const double enph = 1e5;  //KeV (100 MeV) 
  
  const    int Ebin =  50; //KeV (bin of 5 MeV) 
  const    int Tbin =  100; //Time bin in seconds
  static const double de = pow(emax/emin,1.0/Ebin);
  //  static const double dt   = tmax/(Tbin-1);
  const double erg2meV   = 624151.0;


  const double GBM1=10.0;                     // 10 keV 
  const double GBM2=25.0e3;                   // 25 MeV 
  const double LAT1=50.0e3;                   // 50 MeV 
  const double LAT2=3.0e6;                    // 300 GeV 
  const double EGRET1=30.0e3;                 // 30 MeV
  const double EGRET2=100.0e3;                // 100MeV
  const double EGRET3=30.0e6;                  // 30 GeV

  //////////////////////////////////////////////////
  
};

#endif
