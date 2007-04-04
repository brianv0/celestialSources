/*!
 * \class GRBobsengine
 *
 * This class generates a sequance of pulses stored into a std::vector.
 * 
 * \author Nicola Omodei       nicola.omodei@pi.infn.it 
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
  /// Initialize the simulation, getting the burst direction. 
  GRBobsengine(GRBobsParameters *params);
  ~GRBobsengine(){;}
  /*!
    Generate a sequance of pulses (GRBobsPulses) and is store them into a vector.
  */
  double generatePulses(std::vector<GRBobsPulse*> &thePulses, double duration);
  std::vector<GRBobsPulse*> CreatePulsesVector();
  /// Returns the galactic position of the burst in <em>l and <em>b coordinates.
  inline std::pair<double,double> GetDirection(){return m_dir;}
  
 private:
  std::pair<double,double> m_dir;
  GRBobsParameters *m_params;
};

#endif
