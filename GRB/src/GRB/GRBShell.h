/*!
 * \class GRBShell 
 * \brief Describes a shell produced by the blast of the GRB inner engine.
 *
 * Different geometries can be used: jet and issotropic
 *
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
   * \arg The Lorentz factor \c m_gamma is set by GRBengine). 
   * \arg The mass of the shell is computed as: \f$\displaystyle{E\over \Gamma c^2}\f$.
   * .
   * \param E energy of the shell (in erg). In practice, all created shells share
   * the same fraction of the total energy released by the inner engine. 
   */
  GRBShell(double gamma, double mass,
	   double thickness, double radius,double distance = 0.0,std::string type = "jet");
  
  ~GRBShell() { }
  //! This operator compute the resulting shell after the collision between two shells.
  GRBShell operator+(GRBShell Sh2);
  // Accessors to data member values:
  inline double getMass()      {return m_mass;}

  inline double getGamma()     {return m_gamma;}

  inline double getThickness() {return m_thickness;}

  inline double getRadius()    {return m_radius;}  

  inline double getDistance()    {return m_distance;} 

  inline double getBeta()      {return beta(m_gamma);}

  inline double getEint()      {return m_eint;}
  /*! Returns the comoving volume (in \f$cm^3\f$).
   * The comoving volume is defined as 
   *\f$\Large{4\pi\times{thickness}\times{distance}^2}\f$ for isotropic shells,
   *\f$\Large{\pi\times{thickness}\times{radius}^2}\f$ for jet shells.   
   */
  inline double getVolCom()    
    {
      return m_volume;
    }
  
  //Set functions: Should be more protected!!

  inline void setMass(double value)      {m_mass = value;}

  inline void setGamma(double value)     {m_gamma = value;}

  inline void setThickness(double value) {m_thickness = value;}

  inline void setRadius(double value)    {m_radius = value;}

  inline void setDistance(double value)  {m_distance = value;}

  inline void setEint(double value)      {m_eint = value;}
  
  //Higher  level functions:

  /*! 
   * \retval  beta = \f$\sqrt{1-{1\over\Gamma^2}}\f$
   * \param gamma (Lorentz factor of the shell) \f$\Gamma\f$
   */
  double    beta(const double gamma);

  /*!
   * \brief Time evolution of the shell.
   *
   * This method can be used by GRBengine to evolve the shells prior to checking
   * for new shocks.
   * Radius of the shell is moved forward by \f$\beta\times c\times time\f$
   * \param time timestep of evolution
   */
  void      evolve(double time);
    
 private:
  //Data Members:
  double m_gamma;
  double m_mass;
  double m_thickness;
  double m_radius;
  double m_distance;
  double m_volume;
  double m_eint;
};

#endif

