#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include <vector>
#include <string>
#include <sstream>
#include <iostream>

/*! 
 * \namespace TmpCst
 * \brief Namespace containing the constants of the model such as the binning in time and energy
 
 * \author Nicola Omodei       nicola.omodei@pi.infn.it 
 */

namespace TmpCst
{
 /// minimum energy for the computed spectrum (keV)
  const double emin = 10.0;
  /// maximum energy (keV)
  const double emax = 3e8;
  /// default value for the minimum energy of the extracted photons -LAT photons- (keV)
  const double enph = 3e4; 
  //////////////////////////////////////////////////
  /// minimum energy for the computed spectrum (keV)
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
  const double SWIFT1=15.0;
  const double SWIFT2=350.0;      
  //////////////////////////////////////////////////
};

/*!
  \class GRBtemplateSim
  \brief Simulator engine of a GRB source simulated with the GRB phenomenological model.
 
  \author Nicola Omodei       nicola.omodei@pi.infn.it 
 
*/

#ifndef GRBtemplateSIM_H
#define GRBtemplateSIM_H 1


class GRBtemplateSim 
{
 public:

  GRBtemplateSim(std::string InputFileName);
  //! destructor
  ~GRBtemplateSim()
    {
      if(!m_Nv) delete m_Nv;
    }
  /// This method ensures that a unique name is given to the ROOT objects. It is set equal to the pointer address.
  void GetUniqueName(void *ptr, std::string & name);
    
  /*!
   * \brief Starts the GRBtemplate simulation
   *
   * Initialize the simulation
   */
  TH2D* MakeGRB();

  /*! 
    Compute the Flux, as a function of time. It returns a matrix.
  */
  ///    Converts \f$ph/(m^2~s~keV)\f$ into \f$ph/(m^2)\f$, by multiplying for the bin widths.
  TH2D *GRBtemplateSim::Nph(const TH2D *Nv);
  /// Returns the (l,b) galactic position of the burst
  //  inline std::pair<double,double> GRBdir(){return m_GRBengine->GetDirection();}
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
  void GetGBMFlux(std::string GRBname);
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
  void SaveGBMDefinition(std::string GRBname, double ra, double dec, double theta, double phi, double tstart);
  void ComputeEnergyBins();
 private:
  
  /// Gathers all relevant constants for the simulation 
  std::string m_InputFileName;
  double m_tfinal;
  double m_emin;
  double m_emax;
  
  double m_TimeBinWidth;
  int    m_EnergyBins;
  int    m_TimeBins;

  TH2D *m_Nv;
  std::vector<double> NaIEnergyGrid_Vector;
  std::vector<double> BGOEnergyGrid_Vector;
};

#endif


