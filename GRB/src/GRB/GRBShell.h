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
  double GetBeta()  {return m_b;}
  double GetRadius(){return m_r;}
  void   SetRadius(double r){m_r=r;}
  //  double GetTime(){return m_t;}
  //  void   SetTime(double t){m_t=t;}
  double GetGamma(){return m_g;}
  double GetThickness(){return m_dr;}
  double GetMass(){return m_m;} //
  double GetEnergy(){return m_e;}
  double GetVolume();
  double GetComovingVolume(){return GetVolume()*m_g;} // cm^3
  double GetComPartDens();
 private:
  double m_g, m_r,m_dr,m_e, m_m, m_b;
  
};
//////////////////////////////////////////////////
#endif
