/*!
 * \class GRBICompton
 * \brief ICompton emission process.
 * 
 * This class compute the Inverse Compton scattering emission. 
 * This process is represented by \f$\gamma+e\rightarrow\gamma+e\f$
 * The energy of the incident electron (\f$\gamma mc^2\f$) is greater then 
 * the energy of the incident photon (\f$h\nu\f$). The scattered electron give 
 * part of its energy to the photon. The energy of the scattered photon roughly 
 * increases by a factor \f$\gamma^2\f$, unless the Klein Nishina limit is reached.
 * 
 * The easiest way to compute the IC emission is to adopt a power law assumption, 
 * where only a part of the entire population of electrons is used. 
 * In this class either the power law approximation and the 
 * full integration on the electron distribution are used.
 * 
 * \author Nicola Omodei       nicola.omodei@pi.infn.it 
 * \author Johann Cohen-Tanugi johann.cohen@pi.infn.it
 *
 */
#include <vector>
#include"GRBConstants.h"
#include "SpectObj.h"
#include "RadiationProcess.h"
#include "GRBShock.h"

#ifndef GRBICOMPTON_H
#define GRBICOMPTON_H 1

class GRBICompton : virtual public RadiationProcess
{
 public:
  GRBICompton();
  GRBICompton(SpectObj spectrumObj);
  ~GRBICompton(){;}

  /*! Given a Shock compute the spectrum at a certain time
    
    /param *Shock contains all the information to compute the emission
    /param time is the rtime to compute the spectrum
    /param angle is needed to compute the observed spectrum
    /param distance_to_source is the distance of the source (in cm)
   */
  void load(const GRBShock *Shock,
	    const double time  = 0.0,
	    const double angle = 0.0,
	    const double distance_to_source =0.0);

  /*! General Interface for computing the IC spectrum at a certain 
   */
  void load(const double time      = 0.0,
	    const double angle     = 0.0,
	    const double distance_to_source =0.0,
	    const double GAMMAF    = 500.,
	    const double B         = 1.0e+6,
	    const double N0        = 1.0e+50,
	    const double gamma_min = 1.,
	    const double gamma_max = 1.0e+10,
	    const double Vol       = 1.0e+40,
	    const double dr        = 1.0e+10);
  /*! \brief Inverse Compton Function
   *
   This Function represents the probability to obtain a photon photon with 
   energy \e energy from the scattering of an electron of energy
   \f$g_0 m_e c^2\f$ and a photon with energy \f$ e_0 \f$.
   
   \param energy, final energy of the photon after the scattering
   \param e0 initial energy of the photon
   \param g0 Lorentz factor of the incident electron
  
   This function has been derived by Blumenthal and Gould (1970)
   To obtain the Inverse vompton emissivity we have to integrate over
   the distribution of the seed photon and over the electron distribution.
  */
  double InverseComptonFunction(double g0,double e0,double energy);
 private:
};
#endif


