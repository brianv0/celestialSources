/*!
 * \class GRBICompton
 * \brief ICompton emission process.
 * 
 * 
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

  /*! Given a Shock compute the spectrum at a certain /param time.
   *  /param angle is needed to compute the observed spectrum.
   */
  void load(const GRBShock *Shock,
	    const double time  = 0.0,
	    const double angle = 0.0,
	    const double distance_to_source =0.0);

  /*! General Interface for computing the IC spectrum at a certain /param time.
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
   * This Function represents the probability to obtain a photon with energy 
   * \param energy, from the IC scattering between a photon of energy 
   * \e e0 and an electron of energy \f& g0 m_e c^2\f$.  
   * This function has been derived by Blumenthal and Gould (1970)
   * To obtain the Inverse vompton emissivity we have to integrate over
   * the distribution of the seed photon and over the electron distribution.
   */
  double InverseComptonFunction(double g0,double e0,double energy);
 private:
};
#endif


