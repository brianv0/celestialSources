#ifndef GRBSpectrum_H
#define GRBSpectrum_H

#include <vector>
#include <string>
#include <map>
#include <cmath>
#include "flux/ISpectrum.h"
#include "facilities/Observer.h"
#include "flux/GPS.h"
#include "GRBSim.h"
#include "SpectObj.h"

/*!
 * \class GRBSpectrum
 *
 * \brief Main interface to FluxSvc for the physical model
 *  This class inherits from ISpectrum, and allows to pass the result of GRB
 *  physical model simulation to the GLAST simulation framework
 * 
 * \author Nicola Omodei       nicola.omodei@pi.infn.it 
 * \author Johann Cohen-Tanugi johann.cohen@pi.infn.it
 *
 */
class GRBSpectrum : public ISpectrum
{
  friend class GRBmanager;  
 public:
  /*! 
    This initializes the simulation parsing the parameters.
    
    \param params are set in the xml source library in xml directory.
    It contains the name of the file in which the seed is read.
    The seed is necessary to initialize the random number generator.
    Saving the seed in an external file it is possible to reproduce a 
    particular GRB or increase the seed each time a simulation is running.

    An example the xml source declaration for this spectrum should appears:
    \verbatim
    <source name=" GRBspectrum_Gal">
    <spectrum  escale="GeV"> <SpectrumClass name="GRBSpectrum" params="temp.txt"/>
    <!--This will live the direction free-->
    <use_spectrum frame="galaxy"/>
    </spectrum> </source>
    \endverbatim
  */

  GRBSpectrum(const std::string& params);
 
  virtual  ~GRBSpectrum()
    {
      delete m_grbsim;
    }  
  /*! Computes the flux, in \b photons/m^2/s, for a given time
   *  time is actually not used because the spectrum is already updated.
   */
  double flux(double time)const;
  
  /*! \brief Returns the time interval
   *
   * Given \f$t_0\f$, it computes \f$t_1\f$ for which \f$\int_{t_0}^{t_1}1/Rate(t)dt=1\f$
   * and returns the time interval \f$t_1-t_0\f$
   */
  double interval(double time);

  //! Galactic direction 
  inline std::pair<double,double> 
    dir(double){return m_grbsim->GRBdir();} 

  /*! \brief Draws from the current spectrum the energy of a sampled photon. 
   *  \param u uniform random number drawn in the method \c energySrc .  
   */ 
  float operator() (float u ) const ;  
  /*! \brief returns the energy of a sampled photon.
   *
   *  Method called by \c FluxSource::event(). 
   *  It returns the energy of a sampled photon, by calling 
   *  the operator() method. 
   *  \param time    current time. 
   */
  double energy(double time);

  std::string title() const {return "GRBSpectrum";}

  const char * particleName() const {return "gamma";}

  const char * nameOf() const {return "GRBSpectrum";}
  
 private:
  
  //! GRBSim is responsible for the instantiation of the GRB model.
  GRBSim*              m_grbsim;
  
  /*! Binned Spectrum obtained by the GRB model. 
   *  The energy binning is defined in \c GRBConstants.
   *  It is recomputed by \c GRBSim at each time step.
   */ 
  SpectObj m_spectrum;
};
#endif
