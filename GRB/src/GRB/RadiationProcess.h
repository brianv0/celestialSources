/*!
 * \class RadiationProcess
 *
 * \brief This class describes a generic radiative process
 * 
 * Synchrotron Radiation and Inverse Compton scattering inherit from here
 *
 * \author Nicola Omodei       nicola.omodei@pi.infn.it 
 * \author Johann Cohen-Tanugi johann.cohen@pi.infn.it
 *
 */
#ifndef RadiationProcess_H
#define RadiationProcess_H 1

#include "SpectObj.h"
#include"GRBConstants.h"

class RadiationProcess
{

 public:
  RadiationProcess(){;}
  RadiationProcess(SpectObj);
  virtual ~RadiationProcess(){;}
  /*! \brief Represents the power law description for a GRB spectrum.
   *
   * The value it returns is adimensional.
   * \param E is the energy which the flux is returned
   * \param ec is the cooling energy. (depends on the process)
   * \param em is the minimum energy. (depends on the process)
   *
   * Two different regimes are considered :
   * The fast cooling regime appends when \e ec is less then \e em. This means
   * that all the particle can rapidly cool by this process.
   * The slow cooling regime is when \e ec is greater the \e em, 
   * and not all the particles can fast radiate their energy.
   */
  virtual double processFlux(double E, 
			     double ec, 
			     double em);
  
  
  /*! \brief Evolution of the number of electron of given energy.
    This method describes the evolution of the number of electron 
    of energy \f$gi m c^2\f$ in function of the \e ComovingTime.
    The particles have a power law spectral distribution (of index \e p)
    from \e gamma_min to \e gamma_max. The thickness of the region they 
    live is \e dr. Their total number is \e N0.
    A particular cooling time has to be assign at the particular radiation process.
    The cooling time corresponds to the time that a particle takes to emit its energy 
    with this type of radiation process. Different evolution of the number of particle with 
    respect to the time are considered depending on the variable cst::pulse_shape
  */  
  virtual double electronNumber(double gi,
				double gamma_min, 
				double gamma_max,
				double dr,
				double ComovingTime, 
				double CoolingTime,
				double N0);
  /*! \brief Calculates the time shift due to Quantum Gravity
   *
   * The theory of QG preview a dispersion law for the observed photons 
   * that depends on the distance where thay have been produced, 
   * on their observed energy, and on an energy scale that is of the order 
   * of the Plank Energy. In formula:
   *
   * \f$ \delta t = \frac{D E}{c E_{QG}}\f$
   * \param time Time without the effect of QG
   * \param E Observed photon energy (eV)
   * \param distance_to_source Distance of the source, in cm
   * \retval shifted_time  \f$= time -\delta t\f$
   */
  double timeShiftForDispersion(const double time, 
				const double E, 
				const double distance_to_source);

  /*! \brief Calculates the time in the frame where the radiation is emitted. 
   *
   * This method calculates the Lorentz transformation for the time from a 
   * frame that is moving with the radiative region (wich is moving with a
   *  Lorentz factor \e gamma) and our reference frame. 
   * It also consider the time shift due to QG effect.
   
   * \param time Time in the comoving frame of the emitting region
   * \param gamma Lorentz factor of the emitting region
   * \param E Observed photon energy (eV)
   * \param distance_to_source Distance of the source, in cm
   * \retval shifted_time  = gamma * timeShiftForDispersion
   */
  double comovingTime(const double time, 
		      const double gamma, 
		      const double E, 
		      const double distance_to_source);
  
  //! return the Spectrum object.
  inline SpectObj getSpectrumObj() {return m_spectrumObj;}
  inline SpectObj setSpectrumObj(SpectObj spectrum_IN){m_spectrumObj=spectrum_IN;}

  //! Helper function. /param B is in Gauss, Umag is in erg.
  inline virtual double Umag(const double B)
    {
      return pow(B/1.0e+4,2.)/(200.0*M_PI*cst::mu0);
    }
 
 protected:
  SpectObj m_spectrumObj;

};

#endif;
