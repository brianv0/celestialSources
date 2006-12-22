/*!
 *\class  GRBSpectrum
 *
 * Spectrum class for GRB Source Simulation.
 * Inherits from the Spectrum class.
 *
 * \author Nicola Omodei       nicola.omodei@pi.infn.it 
 * \author Johann Cohen-Tanugi johann.cohen@pi.infn.it
 *
 */

#ifndef GRBSpectrum_H
#define GRBSpectrum_H

#include <vector>
#include <string>
#include <map>
#include <cmath>
#include "src/Spectrum.h"
#include "facilities/Observer.h"
#include "src/GPS.h"
#include "CLHEP/Random/RandomEngine.h"
#include "GRBSim.h"

//! Class interfacing the framework with the GRB generation.
class GRBSpectrum : public Spectrum
{
 public:

  //! Constructor: takes a file with some parameters as argument.
  /*! \param params file of parameters 
   */
  GRBSpectrum(const std::string& /*params*/);

  //! Destructor
  ~GRBSpectrum();
  
  //! Computes the flux, in \b photons/m^2/s, for a given time
  double flux(double time)const;

  //! returns rate, for a given time;
  double rate(double time)const;

  /*! \brief Returns the time interval
   *
   * Given \f$t_0\f$, it computes \f$t_1\f$ for which \f$\int_{t_0}^{t_1}1/Rate(t)dt=1\f$
   * and returns the time interval \f$t_1-t_0\f$
   */
  double interval(double time)const;
  // {return 1/rate(time);}

  //! returns the solid angle spanned by the source: set to 1.0 for GRBs.
  double solidAngle() const;
  
  //! Galactic direction 
  std::pair<float,float> dir(float energy) const;
  
  /*! \brief Draws from the current spectrum the energy of a sampled photon. 
   *  \param u uniform random number drawn in the method \c energySrc .  
   */ 
  float operator() (float /*u */ ) const ;
  
  /*! \brief returns the energy of a sampled photon.
   *
   *  Method called by \c FluxSource::event(). 
   *  It returns the energy of a sampled photon, by calling 
   *  the \c operator() method. 
   *  \param engine  random engine for uniform sampling;
   *  \param time    current time. 
   */
  double energySrc(HepRandomEngine*, double /*time*/ );
  
  //! inherited from Spectrum
  inline std::string title() const {return "GRBSpectrum";}
  
  //! inherited from Spectrum
  inline const char * particleName() const {return "gamma";}
  
  //! inherited from Spectrum
  inline  const char * nameOf() const {return "GRBSpectrum";}
  
  
 private:
  
  //! GRBSim is responsible for the instantiation of the GRB model.
  GRBSim*              m_grbsim;
  
  /*! Binned Spectrum obtained by the GRB model. 
   *  The energy binning is defined in \c GRBConstants.
   *  It is recomputed by \c GRBSim at each time step.
   */ 
  std::vector<double>  m_spectrum;
};
#endif
