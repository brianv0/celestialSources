//
//    GRBShell: Class that describes the a Shell
//    Authors: Nicola Omodei & Johann Cohen Tanugi 
//

#include "GRBShell.h"

using namespace cst;

GRBShell::GRBShell(double g, double r, double dr, double e)
{
  m_g  = g;
  m_r  = r;
  m_dr = dr;
  m_e  = e;
  m_m  = e/(g*c2);
  //  m_t  = t;
  m_b=sqrt(1.0 - 1.0/(m_g*m_g));
}

GRBShell::GRBShell(double g, double r, double dr, double e, double m)
{
  m_g  = g;
  m_r  = r;
  m_dr = dr;
  m_e  = e; 
  m_m  = m;
  //  m_t  = t;
  m_b=sqrt(1.0 - 1.0/(m_g*m_g));
}

void GRBShell::Evolve(double t)
{
  m_r  += c * GetBeta() * t;
  //  m_dr += c * GetBeta() * t;
}
//////////////////////////////////////////////////
double GRBShell::GetVolume()
{
  return 4.*cst::pi * m_r * m_r * m_dr;
}  

double GRBShell::GetComPartDens() 
{
  return (m_e * cst::erg2meV)/(m_g*cst::mpc2)/GetComovingVolume(); //1/cm^3
}
