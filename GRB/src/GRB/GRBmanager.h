/*!
  \class  GRBmanager
  
  \brief Spectrum class for many GRBs 
  This class concatenates several GRB one after the other 
  for simulating a series of several GRBs.
  
  \author Nicola Omodei       nicola.omodei@pi.infn.it 
  \author Johann Cohen-Tanugi johann.cohen@pi.infn.it
*/

#ifndef GRBmanager_H
#define GRBmanager_H
#include "GRBConstants.h"

#include "TString.h"

#include <vector>
#include <string>
#include <cmath>
#include "flux/ISpectrum.h"
#include "flux/EventSource.h"
//#include "facilities/Observer.h"
#include "GRBSim.h"
#include "SpectObj/SpectObj.h"

#include "facilities/Util.h"

//class ISpectrum;

class GRBmanager : public ISpectrum
{
  
 public:
  /*! This initializes the simulation parsing the parameters.
    
    \param params are set in the xml source library in xml directory.
    They are: 
    - The time of the first burst
    - The time to wait before the next burst
    - The Minimum photon energy
    
    An example the xml source declaration for this spectrum should appears:
    \verbatim
    <source name=" GRBmanager_Gal">
    <spectrum escale="GeV"> <SpectrumClass name="GRBmanager" params="50 100"/>
    <use_spectrum frame="galaxy"/> 
    </spectrum> </source>
    \endverbatim
  */
  
  GRBmanager(const std::string& params);
  
  virtual  ~GRBmanager();
   
  /*! If a burst is shining it returns the GRBSpectrum::flux method 
   */
  double flux(double time)const;
  /*! \brief Returns the time interval
   *
   * If a burst is shining it returns the GRBSpectrum::interval method.
   * If not it returns the time to whait for the first photon of the next burst.
   */
  double interval(double time);
  
  //! direction, taken from GRBSim
  inline std::pair<double,double> dir(double energy) 
    {
      return m_GalDir;
    } 
  //! calls GRBSpectrum::energySrc
  double energy(double time);

  std::string title() const {return "GRBmanager";} 
  const char * particleName() const {return "gamma";}
  const char * nameOf() const {return "GRBmanager";}
  TString GetGRBname(double time);
  void GenerateGRB();
 
  /*! 
    This method parses the parameter list
    \param input is the string to parse
    \param index if the position of the parameter in the input list. 
    \retval output is the value of the parameter as float number.
  */  
  double parseParamList(std::string input, int index);  
  
 private:
  double m_Rest;
  double m_Frac;
  double m_ra;
  double m_dec;
  double m_l;
  double m_b;
  double m_theta;
  double m_phi;
  
  GRBSim   *m_GRB;
  SpectObj *m_spectrum;
  Parameters *m_par;
  
  const std::string& m_params;
  std::string paramFile;

  double m_timeToWait;
  double m_startTime;
  double m_fluence;
  UInt_t m_GRBnumber;
  double m_enph;
  double m_endTime;

  double m_nextBurst;
  int    m_Nbursts;
  std::pair<double,double> m_GalDir;
  bool m_GenerateGBMOutputs;
};
#endif
