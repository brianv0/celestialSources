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
}

GRBShell::GRBShell(double g, double r, double dr, double e, double m)
{
  m_g  = g;
  m_r  = r;
  m_dr = dr;
  m_e  = e; 
  m_m  = m;
}

void GRBShell::Evolve(double t)
{
  m_r  += c * GetBeta() * t;
  //  m_dr += c * GetBeta() * t;
}
//////////////////////////////////////////////////
