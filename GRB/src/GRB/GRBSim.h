/*!
  \class GRBSim
  \brief Simulator engine of a GRB source.
  
  This class initializes the simulation generating shells with random 
  Lorentz factors and stacks the shells in a vector.
  The evolution with the time is evaluated for all the shells and when two 
  shells have the same radius they collide, and a GRBShock is created. 
  All the shocks are staked in a vector. For each shock is evaluated the 
  observed time at which it occurs.
  This time is used to sort the vector. GRBsim compute the spectrum at a 
  given time and returns it in a vector. 
  This class returns the energy of a photon (chosen from the spectrum), 
  the flux and the rate. 
 
  \author Nicola Omodei       nicola.omodei@pi.infn.it 
  \author Johann Cohen-Tanugi johann.cohen@pi.infn.it
 
*/

#include "GRBShock.h"
#include "GRBConstants.h"
#include "GRBengine.h"
#include "SpectObj.h"
#include <vector>

#ifndef GRBSIM_H
#define GRBSIM_H 1

class GRBSim 
{
 public:
  
  //! The simulation can be initialized by setting the seed of the Random engine generator. If seed is 0, a Random GRB is computed.
  GRBSim(Parameters *params);
  //! destructor
  ~GRBSim(){;}
  
  
  /*!
   * \brief Starts the GRB simulation
   *
   * Initialize the simulation
   */
  TH2D* Fireball();
  
  /*! Compute the Flux, as a function of time. It returns a matrix.
   * \param time is the time in which the spectrum is calculated.
   */
  TH2D *GRBSim::Nph(const TH2D *Nv);
  inline std::pair<float,float> GRBdir(){return m_GRB->GetDirection();}
  inline double Tmax(){return m_tfinal;}
  void SaveNv();
 private:
  
  //! Gathers all relevant constants for the simulation 
  Parameters *m_params;
  GRBengine  *m_GRB;
  double m_tfinal;
  TH2D *m_Nv;
  //std::pair<float,float> m_direction;
  //data member
};

#endif


