/*!
 * \class GRBShock
 *
 * \brief This class implements the shock between 2 shells.
 * 
 * All the parameters that are needed to calculate the emission are 
 * calculated here. 
 * This class calculates the shell which is the result after the shock, the magnetic field and the parameters that determine the distribution of the accelerated electrons.
 * Than it calculates the synchrotron and the inverse Compton emission.   
 *
 * \author Nicola Omodei       nicola.omodei@pi.infn.it 
 * \author Johann Cohen-Tanugi johann.cohen@pi.infn.it
 *
 */

#include <vector>
#include "GRBShell.h"

#ifndef GRBSHOCK_H
#define GRBSHOCK_H 1


class GRBShock 
{
  
 public:

  /*!
   * \brief Constructor: computes the results of the shock between 2 shells.
   *
   * The constructor removes the second shell after having added its content 
   * to the first one. It also initializes several variables.
   * \param Sh1 front shell, the only one remaining after shock  
   * \param Sh2 back shell, deleted after shock.
   * \param time instant of the shock in the local frame
   */
  GRBShock(GRBShell, GRBShell, double time);
  GRBShock(GRBShell);
  //! destructor
  ~GRBShock() { } 
  
  //Accessors:
  //! The time is evaluated in the Shell reference frame
  //  inline double time() {return m_time;}
  //! This is the time seen by GLAST
  inline double tobs() const {return m_tobs;}
  //! Radius of the resulting Shell 
  //  inline double Radius() const {return m_radius;}
  /*! \brief Internal energy. 
   *
   * Is calculated as the difference of the 
   * kinatic energies before and after the collision:
   *\f$E_{int}=m_1c^{2}(\gamma_1-\Gamma_f)+m_2c^{2}(\Gamma_2-\Gamma_f)\f$
   */  
  inline double Eint() const {return m_Eint;}
  /*! \brief The Lorentz factor af the resulting Shell.
   *
   *It is:
   * \f$\Gamma_f\approx\sqrt{\frac{m_1\Gamma_1+m_2\Gamma_2}{m_1/\Gamma_1*m_2/\Gamma_2}}\f$
   */
  inline double getGammaf() const {return m_gf;}
  
  /*! \brief Returns the comoving volume associated to the resulting Shell.
   *
   * This methods calls the VolCom() method of GRBShell
   */
  inline double getVolume() const {return m_volume;}
  inline double getThickness() const {return m_thickness;}

  /*!\brief Is the Total number of barions.
   * 
   * Defined as \f$E_{int}/\gamma m_{p}c^{2}\f$
   * The density of the particles \f$(n_p)\f$ is the ratio 
   * between the total number of and the comoving volume.
   */
  inline double getParticleN() const {return m_partnumber;}

  /*!\brief The Magnetic Field (in Gauss).
   * 
   * The dimensionless parameter \f$\epsilon_B\f$ mesaures the ratio
   * of the magnetic field energy density and the total internal energy.
   * After some simplification:
   * \f$B_{eq}\propto (\epsilon_B)^{1/2}(n_p)^{1/2}*\Gamma_f\f$
   *
   */
  inline double getB() const {return m_Beq;}
  inline double getGammaMin() const {return m_gemin;}
  inline double getGammaMax() const {return m_gemax;}

    //Set functions:

  /*! \brief Set the observer time (in the GLAST frame)
   */
  inline void setTobs(double value) {m_tobs = value;}
  //  inline void setSum(double value)  {m_sum=value;}


  //high level functions:

  //! A printout utility.
  void Write();
  
  //! Return the approximative duration of the shock
  double duration();
  
  //!Returns the flux at time \param time and at \param energy.If \param flagIC if =0 only the syncrothron component is taken into account
  //  std::vector<double> GRBShock::FluxAtT(double time);
  
 private:
  //  double m_time;
  /*! \brief When the time is set also the observed time is calculated. 
   * 
   * The relation between the times in the two different reference frame is:
   * \f$t_{obs}=t-r(t)/c\f$
   */
  double m_tobs;
  double m_mass;
  double m_thickness;
  double m_volume;
  double m_Eint;
  double m_gf;

  double m_partdensity;
  double m_partnumber;
  
  double m_gemin;
  double m_gemax;
  double m_Beq;
  double m_Ue;
  double m_Ub;
  
  /*
    // double m_tsyn;
    double m_riset;
    
    double m_gsh;
    double m_gr;
    double m_gi;
  */
};
#endif
