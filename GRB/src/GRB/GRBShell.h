/*!
 * \class GRBShell 
 * \brief Describes a shell produced by the blast of the GRB inner engine.
 *
 * Different geometries can be used: jet and isotropic
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
  GRBShell(double g, double r, double d, double e);
  GRBShell(double g, double r, double d, double e, double m);
  ~GRBShell(){;}
  
  void Evolve(double dt);
  double GetBeta()  {return sqrt(1.0 - 1.0/(m_g*m_g));}
  double GetRadius(){return m_r;}
  double GetGamma(){return m_g;}
  double GetThickness(){return m_dr;}
  double GetMass(){return m_m;} //
  double GetEnergy(){return m_e;}
  double GetVolume(){return 4./3.*cst::pi*(pow(m_r+m_dr,3.)-pow(m_r,3.));}  // cm^3
  double GetComovingVolume(){return GetVolume()*m_g;} // cm^3
  double GetComPartDens() {return m_m*cst::c2*cst::erg2meV/(cst::mpc2*GetComovingVolume());} //1/cm^3
 private:
  double m_g, m_r,m_dr,m_e, m_m;

};
//////////////////////////////////////////////////
#endif
