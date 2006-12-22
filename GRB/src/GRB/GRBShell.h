/*!
 * \class GRBShell 
 * \brief Describes a spherical shell produced by the blast of the GRB inner engine.
 *
 * The shells will be stacked in a vector.
 * \author Nicola Omodei       nicola.omodei@pi.infn.it 
 * \author Johann Cohen-Tanugi johann.cohen@pi.infn.it
 *
 */

#include "GRBConstants.h"

#ifndef GRBSHELL_H
#define GRBSHELL_H 1

class GRBShell 
{  
 public:

  /*!
   * \brief Creation of a new shell.
   *
   * \arg The Lorentz factor \c m_gamma is randomly drawn (see \c generateGamma). 
   * \arg The mass of the shell is computed as: \f$\displaystyle{E\over \Gamma c^2}\f$.
   * .
   * \param E energy of the shell (in erg). In practice, all created shells share
   * the same fraction of the total energy released by the inner engine. 
   */
  //  GRBShell(double /*energy*/);
  GRBShell();

  ~GRBShell() { }
  
  // Accessors to data member values:
  inline double Mass()      {return m_mass;}
  inline double Gamma()     {return m_gamma;}
  inline double Thickness() {return m_thickness;}
  inline double Radius()    {return m_radius;}

  /*! \brief computes and returns the comoving volume (in \f$cm^3\f$).
   *
   * The comoving volume is defined as 
   *\f$\Large{4\pi\times{thickness}\times{radius}^2\times\Gamma}\f$
   */
  inline double VolCom()    {return 4.0*cst::pi*m_thickness*m_radius*m_radius*m_gamma;}
  
  //Set functions: Should be more protected!!
  inline void setMass(double value)      {m_mass = value;}
  inline void setGamma(double value)     {m_gamma = value;}
  inline void setThickness(double value) {m_thickness = value;}
  inline void setRadius(double value)    {m_radius = value;}
 
  //Higher  level functions:


  /*!
   * \brief generate a random Lorentz factor.
   *
   * For a uniform random number \c u, the method returns \f$\Gamma_0+ud\Gamma\f$.
   *\param gamma0 \f$\Gamma_0\f$
   *\param dgamma \f$d\Gamma\f$
   */
  double    generateGamma(double gamma0, double dgamma);

  /*! 
   * \retval  beta = \f$\sqrt{1-{1\over\Gamma^2}}\f$
   * \param gamma (Lorentz factor of the shell) \f$\Gamma\f$
   */
  double    beta(const double gamma);

  /*!
   * \brief Time evolution of the shell.
   *
   * This method is used in GRBSim to evolve the shells prior to checking
   * for new shocks.
   * Radius of the shell is moved forward by \f$\beta\times c\times time\f$
   * \param time timestep of evolution
   */
  void      evolve(double time);
  
  
 private:
  //Data Members:
  double m_mass;
  double m_gamma;
  double m_thickness;
  double m_radius;
};

#endif

