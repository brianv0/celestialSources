//
//    GRBShell: Class that describes the a Shell
//    Authors: Nicola Omodei & Johann Cohen Tanugi 
//

#include "GRBShell.h"

using namespace cst;

GRBShell::GRBShell(double g, double r, double dr, double e, double t)
{
  m_g  = g;
  m_r  = r;
  m_dr = dr;
  m_e  = e;
  m_m  = e/(g*c2);
  m_t  = t;
  m_b=sqrt(1.0 - 1.0/(m_g*m_g));
}

GRBShell::GRBShell(double g, double r, double dr, double e, double m, double t)
{
  m_g  = g;
  m_r  = r;
  m_dr = dr;
  m_e  = e; 
  m_m  = m;
  m_t  = t;
  m_b=sqrt(1.0 - 1.0/(m_g*m_g));
}

void GRBShell::Evolve(double t)
{
  m_r  += c * GetBeta() * t;
  //  m_dr += c * GetBeta() * t;
}
//////////////////////////////////////////////////
double GRBShell::GetVolume(double time)
{
  double r = GetRadius(time);
  return 4.*cst::pi*r*r*m_dr;
}  

double GRBShell::GetComPartDens(double time) 
{
  return (m_e * cst::erg2meV)/(m_g*cst::mpc2)/GetComovingVolume(time); //1/cm^3
}
