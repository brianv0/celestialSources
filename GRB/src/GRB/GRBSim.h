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
#include "GRBSynchrotron.h"
#include "GRBICompton.h"
#include <vector>

#ifndef GRBSIM_H
#define GRBSIM_H 1

class GRBSim 
{
 public:
  
  //! The simulation can be initialized by setting the seed of the Random engine generator. If seed is 0, a Random GRB is computed.
  GRBSim(long seed=0); 
  GRBSim(const std::string& params);
  //! destructor
  ~GRBSim();
  
  //! Gathers all relevant constants for the simulation 
  GRBConstants *myParam;
  
  /*!
   * \brief Starts the GRB simulation
   *
   * Initialize the simulation
   */
  void MakeGRB(double time_offset = 0.);
  
  /*! Compute the Flux, at given time. It returns a vector.
   * \param time is the time in which the spectrum is calculated.
   */
  SpectObj ComputeFlux(const double time);
  //! Sets the private data member m_spectobj to \e myspectrum.
  inline void setSpectrum(SpectObj myspectrum) 
    {m_spectobj=myspectrum;}

  inline SpectObj getSpectrum() 
    {return m_spectobj;}
  //! The direction of the GRB is chosen randomly in the sky.
  inline std::pair<float,float> GRBdir()       {return m_direction;}
  /*! Is the vector that contains the energy bin,in \f$eV\f$, 
    in wich the flux is evaluated.
   */
  //  inline std::vector<double>    Energy()       {return m_energy;}
  //! Is the energy bin size in \f$eV\f$ (the energy is in log scale!)
  // inline std::vector<double>    DeltaE()       {return m_de;}
  
  //! Corresponds to the time (in the GLAST frame) in which the burst ends
  inline double                 Tmax()         {return m_duration;}
  
  /*! Is the Area of sphere having as radius 
    the distance of the source (in \f$m^2\f$)
   */ 
  inline double                 Area()         {return m_area;}
  
  /*! Return the minimum energy for the photon drawn */
  inline double EnergyPh(){return m_enph;}
  /*! Parse the parameter list. */
  long parseParamList(std::string input, int index);
 private:
  //data member
  std::vector<GRBShell> theShells;
  std::vector<GRBShock> theShocks;
  std::pair<float,float> m_direction;
  
  double m_duration;
  double m_ftot;
  double m_phtot;
  double m_distance;
  double m_area;
  double m_DeadTime;
  double m_enph;
  double m_jetangle;
  long m_seed;
  //  HepRandomEngine *m_engine;
  GRBSynchrotron m_synchrotron;
  GRBICompton m_icompton;
  SpectObj m_spectobj;
};

#endif


