/*!
 * \class GRBShock
 *
 * \brief This class implements the shock between 2 shells.
 * 
 * All the parameters that are needed to calculate the emission are 
 * calculated here. 
 * This class calculates the shell which is the result after the shock,
 * the magnetic field and the parameters that determine the distribution of the accelerated electrons.
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
   * \brief Constructor: computes the shocks physics in the Shocked region.
   *
   * The GRBShell represents the shocked region. 
   * The constructor calculate the magnetic field and the particle acceleration in this region.
   */
  GRBShock(GRBShell);

  ~GRBShock() { } 
  
  //Accessors:
  //! This is the time seen by GLAST
  inline double tobs() const {return m_tobs;}

  /*! \brief Internal energy. 
   *
   * Is the internal energy in the shocked material. 
   * Part of it is converted into magnetic field and 
   * part of it needs to accelerate particles. 
   */  

  inline double Eint() const {return m_Eint;}
  /*! \brief The Lorentz factor af the shocked material.
   *
   *It is needed to compute the energy transformation.
   */
  inline double getGammaf() const {return m_gf;}
  
  inline double getVolume() const {return m_volume;}
  inline double getThickness() const {return m_thickness;}

  inline double getParticleN() const {return m_partnumber;}

  inline double getB() const {return m_Beq;}
  inline double getGammaMin() const {return m_gemin;}
  inline double getGammaMax() const {return m_gemax;}

  //Set functions:
  inline void setTobs(double value) {m_tobs = value;}
  //high level functions:

  //! A printout utility.
  void Write();

  //! Return the approximative duration of the shock
  double duration();
    
 private:
  
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
};
#endif
