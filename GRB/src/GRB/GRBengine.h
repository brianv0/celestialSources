/*!
 * \class GRBengine
 *
 * \brief This class permits to generate several GRB.
 *
 * This class sets the parameters of the model, and create a vector of GRBShock.
 * 
 * \author Nicola Omodei       nicola.omodei@pi.infn.it 
 * \author Johann Cohen-Tanugi johann.cohen@pi.infn.it
 *
 */

#ifndef GRBENGINE_H
#define GRBENGINE_H 1

#include "GRBConstants.h"
#include "GRBShock.h"

class GRBengine
{
  
 public:
  /// Constructor. Get the Parameters class as input configuration.
  GRBengine(Parameters *params);
  ///  Default destructor
  ~GRBengine(){;}
  /*!
    Create a vector of GRBShock. The basic idea is:
    - First it gets the parameters from the class Parameters.
    - Start a loop checking on the duration of the burst
    For a two given shells, one slower (front) and one faster (back), a GRBShok is created and a vector is filled.
    
  */
  std::vector<GRBShock*> CreateShocksVector();
  
  /*
    Returns the luminosity distance. The redshift is fixed to 1 since this method 
    is used only for computing the Quantum Gravity effect.
 */
  double GetDistance();
  /// Returns the galastic (l,b) position of the burst
  inline std::pair<double,double> GetDirection(){return m_dir;}
 private:
  std::pair<double,double> m_dir;
  Parameters *m_params;
};

#endif
