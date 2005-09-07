/*!
 * \class GRBobsengine
 *
 * \brief This class permits to generate different kind of GRBs.
 * 
 * \author Nicola Omodei       nicola.omodei@pi.infn.it 
 * \author Johann Cohen-Tanugi johann.cohen@pi.infn.it
 *
 */

#ifndef GRBobsENGINE_H
#define GRBobsENGINE_H 1

#include "GRBobsConstants.h"
#include "GRBobsPulse.h"
#include <vector>
#include "TRandom.h"

class GRBobsengine
{
  
 public:
  GRBobsengine(GRBobsParameters *params);
  ~GRBobsengine()
    {
      //      delete m_params;
    }
    
  std::vector<GRBobsPulse*> CreatePulsesVector();
  
  /*
    inline double GetRiseTime()    {return m_rise;}
    inline double GetDecayTime()   {return m_decay;}
    inline double GetPulseHeight() {return m_pulseHeight;}
  */
  inline std::pair<double,double> GetDirection(){return m_dir;}

 private:
  std::pair<double,double> m_dir;
  GRBobsParameters *m_params;
  //  double m_rise,m_decay,m_pulseHeight;  
};

#endif
