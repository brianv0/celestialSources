/*!
  \class  GRBtemplateManager
  
  \brief Spectrum class for many GRBs 
  This class concatenates several GRB one after the other 
  for simulating a series of several GRBs.
  
  \author Nicola Omodei       nicola.omodei@pi.infn.it 
*/

#ifndef GRBtemplateManager_H
#define GRBtemplateManager_H

#include <map>
#include <cmath>
#include "TRandom.h"
#include "flux/Spectrum.h"
#include "flux/EventSource.h"
#include "GRBtemplateSim.h"
#include "SpectObj/SpectObj.h"

#include "facilities/Util.h"

class GRBtemplateManager : public Spectrum
{
  
 public:
  /*! This initializes the simulation parsing the parameters.
    
    \param params are set in the xml source library in xml directory.
    They are: 
    - The time of the first burst
    - The time to wait before the next burst
    
    An example the xml source declaration for this spectrum should appears:
    \verbatim
    <source name=" GRBtemplateManager_Gal">
    <spectrum escale="GeV"> <SpectrumClass name="GRBtemplateManager" params="50 100"/>
    <use_spectrum frame="galaxy"/> 
    </spectrum> </source>
    \endverbatim
  */
  
  GRBtemplateManager(const std::string& params);
  virtual  ~GRBtemplateManager();
   
  /*! If a burst is shining it returns the flux method 
   */
  double flux(double time)const;
  /*! \brief Returns the time interval
   *
   * If a burst is shining it returns the interval method.
   * If not it returns the time to whait for the first photon of the next burst.
   */
  double interval(double time);
  
  //! direction, taken from GRBtemplateSim
  inline std::pair<double,double> dir(double) 
    {
      //      (void)(energy);//energy is not used here
      return m_GalDir;
    } 
  
  double energy(double time);
  
  std::string title() const {return "GRBtemplateManager";} 
  const char * particleName() const {return "gamma";}
  const char * nameOf() const {return "GRBtemplateManager";}
  std::string GetGRBname();
  void GenerateGRB();  
  void DeleteGRB();  

  /*! 
    This method parses the parameter list
    \param input is the string to parse
    \param index if the position of the parameter in the input list. 
    \retval output is the value of the parameter as float number.
  */  
  std::string parseParamList(std::string input, unsigned int index);  
  
 private:
  double m_Rest;
  double m_Frac;
  double m_ra;
  double m_dec;
  double m_theta,m_theta_fow;
  double m_phi;


  SpectObj    *m_spectrum;
  GRBtemplateSim   *m_GRB;
  TRandom *m_rnd;

  const std::string& m_params;
  std::string paramFile;
  std::pair<double,double> m_GalDir;
  bool m_grbGenerated;
  bool m_grbdeleted;
  bool m_grbocculted;
  bool m_GenerateGBMOutputs;
  
  std::string m_InputFileName;
  double m_l;
  double m_b;
  //  double m_fluence;
  long   m_GRBnumber;
  double m_MinPhotonEnergy;
  double m_startTime;
  double m_endTime;
};
#endif
