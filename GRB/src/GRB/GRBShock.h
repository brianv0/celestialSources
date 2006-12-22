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
   */
  GRBShock(GRBShell*, GRBShell*);

  //! destructor
  ~GRBShock() { }
  
  //Accessors:
  //! The time is evaluated in the Shell reference frame
  inline double time() {return m_time;}
  //! This is the time seen by GLAST
  inline double tobs() {return m_tobs;}
  //! Radius of the resulting Shell 
  inline double Radius() {return m_radius;}
  /*! \brief Internal energy. 
   *
   * Is calculated as the difference of the 
   * kinatic energies before and after the collision:
   *\f$E_{int}=m_1c^{2}(\gamma_1-\Gamma_f)+m_2c^{2}(\Gamma_2-\Gamma_f)\f$
   */  
  inline double Eint() {return m_Eint;}
  /*! \brief The Lorentz factor af the resulting Shell.
   *
   *It is:
   * \f$\Gamma_f\approx\sqrt{\frac{m_1\Gamma_1+m_2\Gamma_2}{m_1/\Gamma_1*m_2/\Gamma_2}}\f$
   */
  inline double gf() {return m_gf;}
  
  /*! \brief Returns the comoving volume associated to the resulting Shell.
   *
   * This methods calls the VolCom() method of GRBShell
   */
  inline double VolCom() {return m_VolCom;}

  /*!\brief Is the Total number of barions.
   * 
   * Defined as \f$E_{int}/\gamma m_{p}c^{2}\f$
   * The density of the particles \f$(n_p)\f$ is the ratio 
   * between the total number of and the comoving volume.
   */
  inline double npart() {return m_npart;}

  /*!\brief The Magnetic Field (in Gauss).
   * 
   * The dimensionless parameter \f$\epsilon_B\f$ mesaures the ratio
   * of the magnetic field energy density and the total internal energy.
   * After some simplification:
   * \f$B_{eq}\propto (\epsilon_B)^{1/2}(n_p)^{1/2}*\Gamma_f\f$
   *
   */
  inline double Beq() {return m_Beq;}
  
  /*!\brief Normalization to the correct emitted flux
   *
   * This metod is used only for normalize the spectrum, and the value is 
   * set by GRBSim.
   */
  inline double Sum() {return m_sum;}

  //Set functions:

  /*! \brief When the time is set also the observed time is calculated. 
   * 
   * The relation between the times in the two different reference frame is:
   * \f$t_{obs}=t-r(t)/c\f$
   */
  inline void setTime(double value) {m_time = value;m_tobs=m_time-m_radius/cst::c;}
  /*! \brief Set the observer time (in the GLAST frame)
   */
  inline void setTobs(double value) {m_tobs = value;}
  inline void setSum(double value)  {m_sum=value;}


  //high level functions:

  //! A printout utility.
  void Write();
  /*!\brief Calculates the synchrotron critical energy (in eV)
   * for an electron.
   *
   * \f$E_{syn}\propto(\Gamma_f B_{eq})(\gamma_e)^2\f$ 
   * \param gammae is the Lorentz factor of the electron \f$\gamma_e\f$ 
    */
  double Esyn(double gammae);
  /*!\brief Calculates the Inverse Compton emission energy (in eV)
   *
   * If \f$h\nu_{in}\f$ is the energy of the incoming photon, and \f$\gamma_e\f$
   * is the Lorentz factor of the incoming electron, the energy 
   * of the resulting photons will be:
   * \f$h\nu_{fin}\approx h\nu_{in}\gamma_e\f$
   */
  double Eic(double gamme);

  /*!
   * Not yet implemented!
   */
  double OpticalDepht();
  
  /*! \brief The synchrotron cooling time.
   *
   * It depends on the energy  of the observed photon.
   * \f$t_{syn}\propto B_{eq}^{-5/2}\sqrt{\Gamma_f/phenergy}\f$
   * \param phenergy is the energy of the observed photon
   */ 
  double tsyn(double phenergy);

  /*! \brief Fast Rise Exponential Decay.
   *
   * This method calculates the temporal behaviour of the spectrum.
   * The temporal profile is well approximated by
   * a double exponential law. In which 
   * the rise time depends on the geometrical size of the shell, 
   * while the decay time depends on the cooling time of the 
   * processes (Syn, IC, ...).
   */
  double fred(double ee, double tt);
  
  //! Calculates the Synchrotron component of the spectrum. 
  double Fsyn(double ee, double tt);
  //! Calculate the Inverse compton part of the spectrum.
  double Fic(double ee, double tt);
  
 private:
  GRBShell* Sh1;
  GRBShell* Sh2;
  double m_time;
  double m_tobs;
  double m_radius;
  double m_mass;
  double m_Eint;
  double m_gsh;
  double m_VolCom;
  double m_npart;
  
  double m_gemin;
  double m_gemax;
  double m_gecoo;
  
  double m_tsyn;
  double m_riset;
  double m_sum;
  
  double m_gf;
  double m_gr;
  double m_gi;
  double m_Beq;
  
};
#endif
