/*!
 * \class GRBShell 
 * \brief Describes a shell produced by the blast of the GRB inner engine.
 *
 * \author Nicola Omodei       nicola.omodei@pi.infn.it 
 * \author Johann Cohen-Tanugi johann.cohen@pi.infn.it
 *
 */

#ifndef GRBSHELL_HH
#define GRBSHELL_HH 1
#include "GRBConstants.h"
//////////////////////////////////////////////////
class GRBShell
{
 public:
  /*!
    Constructor
    \param g: lorentz factor of the shell (\f$\Gamma\f$)
    \param r: radius of the shell, in cm from the central engine.
    \param d: thickness of the shell (in cm)
    \param e: energy pf the shell (in erg)
    the mass of the shell is then computed using the energy-mass relation 
  */
  GRBShell(double g, double r, double d, double e);
  /*!
    Constructor
    \param g: lorentz factor of the shell (\f$\Gamma\f$)
    \param r: radius of the shell, in cm from the central engine.
    \param d: thickness of the shell (in cm)
    \param e: energy pf the shell (in erg)
    \param m: mass of the shell (in grams)
  */
  GRBShell(double g, double r, double d, double e, double m);
  /// destructor
  ~GRBShell(){;}
  /// Evolves the radius of the shell (\f$ \dot{r}=\beta~c\f$)
  void Evolve(double dt);
  double GetBeta()  {return m_b;}
  double GetRadius(){return m_r;}
  void   SetRadius(double r){m_r=r;}
  double GetGamma(){return m_g;}
  double GetThickness(){return m_dr;}
  double GetMass(){return m_m;} //
  double GetEnergy(){return m_e;}
  /*!
    Returns the Volume:
    \f$ V= 4\pi r^2 d \f$
  */
  double GetVolume();
  /*!
    Retuns the comoving volume \f$V_{com}~=~V~\Gamma \f$ 
  */
  double GetComovingVolume(){return GetVolume()*m_g;} // cm^3
  /*!
    Returns the comoving particles dentity: 
    \f$ \frac{e}{\Gamma m_p c^2}\frac{1}{V_{com}}\f$
  */
  double GetComPartDens();
 private:
  double m_g, m_r,m_dr,m_e, m_m, m_b;
  
};
//////////////////////////////////////////////////
#endif
