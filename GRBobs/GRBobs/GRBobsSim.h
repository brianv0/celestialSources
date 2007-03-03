/*!
  \class GRBobsSim
  \brief Simulator engine of a GRB source simulated with the GRB phenomenological model.
 
  \author Nicola Omodei       nicola.omodei@pi.infn.it 
 
*/

#include "GRBobsConstants.h"
#include "GRBobsPulse.h"
#include "GRBobsengine.h"
#include "TString.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#ifndef GRBobsSIM_H
#define GRBobsSIM_H 1

class GRBobsSim 
{
 public:

  GRBobsSim(GRBobsParameters *params);
  //! destructor
  ~GRBobsSim()
    {
      delete m_GRBengine;
      delete m_Nv;
    }
  /// This method ensures that a unique name is given to the ROOT objects. It is set equal to the pointer address.
  void GetUniqueName(const void *ptr, std::string & name);
    
  /*!
   * \brief Starts the GRBobs simulation
   *
   * Initialize the simulation
   */
  TH2D* MakeGRB();

  /*! 
    Compute the Flux, as a function of time. It returns a matrix.
  */
  ///    Converts \f$ph/(m^2~s~keV)\f$ into \f$ph/(m^2)\f$, by multiplying for the bin widths.
  TH2D *GRBobsSim::Nph(const TH2D *Nv);
  /// Returns the (l,b) galactic position of the burst
  inline std::pair<double,double> GRBdir(){return m_GRBengine->GetDirection();}
  /// Rerturn the duration of the burst
  inline double Tmax(){return m_tfinal;} 
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
  
  /// Gathers all relevant constants for the simulation 
  GRBobsParameters *m_params;
  GRBobsengine  *m_GRBengine;
  double m_tfinal;
  double m_fluence;
  int m_tbin;
  TH2D *m_Nv;
};

#endif


