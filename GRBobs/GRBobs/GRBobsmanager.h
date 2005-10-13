/*!
  \class  GRBobsmanager
  
  \brief Spectrum class for many GRBs 
  This class concatenates several GRB one after the other 
  for simulating a series of several GRBs.
  
  \author Nicola Omodei       nicola.omodei@pi.infn.it 
  \author Johann Cohen-Tanugi johann.cohen@pi.infn.it
*/

#ifndef GRBobsmanager_H
#define GRBobsmanager_H
#include "GRBobsConstants.h"

#include "TString.h"

#include <vector>
#include <string>
#include <map>
#include <cmath>
#include "flux/ISpectrum.h"
#include "flux/EventSource.h"
#include "GRBobsSim.h"
#include "SpectralComponent.h"

#include "SpectObj/SpectObj.h"

#include "facilities/Util.h"

//class ISpectrum;

class GRBobsmanager : public ISpectrum
{
  
 public:
  /*! This initializes the simulation parsing the parameters.
    
    \param params are set in the xml source library in xml directory.
    They are: 
    - The time of the first burst
    - The time to wait before the next burst
    
    An example the xml source declaration for this spectrum should appears:
    \verbatim
    <source name=" GRBobsmanager_Gal">
    <spectrum escale="GeV"> <SpectrumClass name="GRBobsmanager" params="50 100"/>
    <use_spectrum frame="galaxy"/> 
    </spectrum> </source>
    \endverbatim
  */
  
  GRBobsmanager(const std::string& params);
  virtual  ~GRBobsmanager();
   
  /*! If a burst is shining it returns the flux method 
   */
  double flux(double time)const;
  /*! \brief Returns the time interval
   *
   * If a burst is shining it returns the interval method.
   * If not it returns the time to whait for the first photon of the next burst.
   */
  double interval(double time);
  
  //! direction, taken from GRBobsSim
  inline std::pair<double,double> dir(double) 
    {
      (void)(energy);//energy is not used here
      return m_GalDir;
    } 
  
  double energy(double time);
  
  std::string title() const {return "GRBobsmanager";} 
  const char * particleName() const {return "gamma";}
  const char * nameOf() const {return "GRBobsmanager";}
  TString GetGRBname();
  void GenerateGRB();  
  void DeleteGRB();  
  /*! 
    This method parses the parameter list
    \param input is the string to parse
    \param index if the position of the parameter in the input list. 
    \retval output is the value of the parameter as float number.
  */  
  double parseParamList(std::string input, unsigned int index);  
  
 private:
  double m_Rest;
  double m_Frac;
  double m_ra;
  double m_dec;
  double m_theta;
  double m_phi;

  SpectObj    *m_spectrum;
  SpectObj    *m_spectrum1;

  SpectralComponent *PromptEmission;
  SpectralComponent *AfterGlowEmission;
  
  GRBobsSim   *m_GRB;
  GRBobsParameters  *m_par;

  const std::string& m_params;
  std::string paramFile;
  std::pair<double,double> m_GalDir;
  bool m_grbGenerated;
  bool m_grbdeleted;
  bool m_GenerateGBMOutputs;

  double m_l;
  double m_b;
  double m_duration;
  double m_fluence;
  double m_z;
  long   m_GRBnumber;
  double m_alpha;
  double m_beta;
  double m_MinPhotonEnergy;
  double m_Energy_CO;
  double m_LATphotons;
  double m_EC_delay;
  double m_EC_duration;
  double m_startTime;
  double m_endTime;
  double m_startTime_EC;
  double m_endTime_EC;


};
#endif
