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

#ifndef GRBSYNCHROTRON_H
#define GRBSYNCHROTRON_H 1

class GRBSynchrotron
{
 public:
  GRBSynchrotron::GRBSynchrotron(){;}
  GRBSynchrotron::GRBSynchrotron(SpectObj spectrumObj);
  
  ~GRBSynchrotron(){;}
  /*! Given a Shock compute the spectrum at a certain /param time.
   *  /param angle is needed to compute the observed spectrum.
   */
  void load(const GRBShock *Shock,
	    const double time  = 0.0,
	    const double angle = 0.0);
  /*! General Interface for computing the synchrotron spectrum at a certain /param time.
   */
  void load(const double time      = 0.0,
	    const double angle     = 0.0,
	    const double GAMMAF    = 500.,
	    const double B         = 1.0e+6,
	    const double N0        = 1.0e+50,
	    const double gamma_min = 1.,
	    const double gamma_max = 1.0e+10,
	    const double Vol       = 1.0e+40,
	    const double dr        = 1.0e+10); 
  //! return the Spectrum object after computing the synchrotron emission.
  inline SpectObj getSpectrumObj() {return m_spectrumObj;}
  //! Helper function. /param B is in Gauss, Umag is in erg.
  inline double Umag(const double B)
    {
      return pow(B/1.0e+4,2.)/(200.0*M_PI*cst::mu0);
    }
  
 private:
  SpectObj m_spectrumObj;
  std::vector<double>    x;
  std::vector<double>    g;
  std::vector<double>   dx;
  std::vector<double> fsyn;
};
#endif


