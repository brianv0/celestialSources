/*!
 * \class GRBSynchrotron
 * \brief Synchrotron emission process.
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
#include "GRBShock.h"

#include "RadiationProcess.h"

#ifndef GRBSYNCHROTRON_H
#define GRBSYNCHROTRON_H 1

class GRBSynchrotron : virtual public RadiationProcess
{
 public:
  GRBSynchrotron();
  GRBSynchrotron(SpectObj spectrumObj);
  
  ~GRBSynchrotron(){;}
  /*! Given a Shock compute the spectrum at a certain /param time.
   *  /param angle is needed to compute the observed spectrum.
   */
  void load(const GRBShock *Shock,
	    const double time  = 0.0,
	    const double angle = 0.0,
	    const double distance_to_source = 0.0);
  /*! General Interface for computing the synchrotron spectrum at a certain /param time.
   */
  void load(const double time      = 0.0,
	    const double angle     = 0.0,
	    const double distance_to_source = 0.0,
	    const double GAMMAF    = 500.,
	    const double B         = 1.0e+6,
	    const double N0        = 1.0e+50,
	    const double gamma_min = 1.,
	    const double gamma_max = 1.0e+10,
	    const double Vol       = 1.0e+40,
	    const double dr        = 1.0e+10); 

  /*! \brief Synchrotron Function
   *
   * This Function represents the approximated version of the Synchrotron 
   * function (eq. [6.31c] in Rybiki & Lightman 1979) and is the normalized 
   * emission of an electron with gyrofrequency \e esyn in function of the 
   * photon \e energy.
   */
  double SynchrotronFunction(double esyn, double energy);
 private:
};
#endif


