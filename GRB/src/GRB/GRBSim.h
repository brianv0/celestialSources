/*!
  \class GRBSim
  \brief Simulator engine of a GRB source.
  
  This class manages the GRB simulation.
  From the Parameters class (which contains the configuration): 
  -# It creates a GRBengine
  -# Gets the BATSE fluence (needed for flux normalization)
  -# Compute the Fireball spectrum, by adding the contribution of all the shocks (GRBShock) created by GRBengine class.
  (The spectrum is in units of [\f$ph/(m^2~s~keV)\f$] and it is computed as a function of time (sec) and energy (keV).
  -# Finally it normalizes the spectrum at BATSE energy: the observed BASTE fluence distribution is obtained.
  
  The GBM output is stored in a definition file, and in a GBM file, which contains the results of a "Band fit" every 16 ms.
   
  \author Nicola Omodei       nicola.omodei@pi.infn.it 
  \author Johann Cohen-Tanugi johann.cohen@pi.infn.it
 
*/
#include "TString.h"
#include "GRBShock.h"
#include "GRBConstants.h"
#include "GRBengine.h"
#include "TH2D.h"
#include <vector>
#include <string>
#include <sstream>
#include <iostream>

#ifndef GRBSIM_H
#define GRBSIM_H 1

class GRBSim 
{
 public:
  
  /// Constructor: the Parameters class is passed as input.
  GRBSim(Parameters *params);
  /// Destructor
  ~GRBSim()
    {
      delete m_GRBengine;
      delete m_Nv;
    }
  /// This method ensures that a unique name is given to the ROOT objects. It is set equal to the pointer address.
  void GetUniqueName(const void *ptr, std::string & name);  
  /// Compute the TH2 object wich stores the flux as a function of energy and time.
  TH2D* Fireball();
  ///    Converts \f$ph/(m^2~s~keV)\f$ into \f$ph/(m^2)\f$, by multiplying for the bin widths.
  TH2D *GRBSim::Nph(const TH2D *Nv);
  /// Returns the (l,b) galactic position of the burst
  inline std::pair<double,double> GRBdir(){return m_GRBengine->GetDirection();}
  inline double GetDistance(){return m_GRBengine->GetDistance();}
  /// Rerturn the duration of the burst
  inline double Tmax(){return m_tfinal;}
  /// Return the fluence of the burst \f$erg/cm^2\f$ in the BATSE energy range (20 keV - 1 MeV);
  inline double GetFluence(){return m_fluence;} 
  /// The GRB number is uses for recognize bursts. It is usually set equal to the GRB starting date.
  inline UInt_t GetGRBNumber(){return m_params->GetGRBNumber();} //erg/cm^2
  /// This methods save the TH2 histogram in a root file. It is used by GRBROOTTest.cxx
  void SaveNv();
  /*!
    \brief This methods performs a series of fits using the well known Band function. (Band et al.(1993) ApJ.,413:281-292)
    
    The flux is fitted every 16 ms and the four parameters of the Band function are saved in a txt file. 
    The name of the txt file is chosen in agreement with the GRB name (See GRBmanager).
    \param GRBname is the name of the GRB created in GRBmanager. It is used for naming the GBM output file, and it is usually computed with the dating convention.
  */
  void GetGBMFlux(TString GRBname);
  /*!
    \brief This methods saves the definition file for GBM simulator.
    
    The standard file format is:
    \verbatim
    std::ofstream os(name,std::ios::out);
    os<<"BURST DEFINITION FILE"<<std::endl;
    os<<"Burst Name"<<std::endl;
    os<<GRBname<<std::endl;
    os<<"RA,DEC (deg):"<<std::endl;
    os<<ra<<" "<<dec<<std::endl;
    os<<"S/C azimuth, elevation (deg):"<<std::endl;
    os<<phi<<" "<<theta<<std::endl;
    os<<"Trigger Time (s):"<<std::endl;
    os<<tstart<<std::endl;
    os.close();
    \endverbatim
    
    Where:
    
    \param GRBname is the string containing the name of the GRB (year,month,day,frac_of_day)
    \param ra is the right ascension of the burst
    \param dec is the declination of the burst
    \param theta is the space craft azimuth (in degree).
    \param phi is the elevation (deg) from LAT horizon (zenith -> phi=90)
    \param tstart is the GRB starting time (in second, since the starting time of the simulation).
  */
 
  void SaveGBMDefinition(TString GRBname, double ra, double dec, double theta, double phi, double tstart);

 private:
  
  Parameters *m_params;
  GRBengine  *m_GRBengine;
  double m_tfinal;
  double m_fluence;
  int Tbin;
  TH2D *m_Nv;
};

#endif


