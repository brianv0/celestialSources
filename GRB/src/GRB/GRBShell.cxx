//
//    GRBShell: Class that describes the a Shell
//    Authors: Nicola Omodei & Johann Cohen Tanugi 
//

#include "GRBShell.h"
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Random/RandExponential.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RanluxEngine.h"

using namespace std;

//GRBShell::GRBShell(double E) 
GRBShell::GRBShell() 
{ 
  //  cout<<cst::gamma0<<endl;
  //m_gamma = generateGamma(cst::gamma0,cst::dgamma);
  //m_mass = E/(m_gamma*cst::c2);
}


double GRBShell::generateGamma(double gamma0,double dgamma) 
{
  HepRandom::setTheEngine(new RanluxEngine);
  double gamma = gamma0 + (double (RandFlat::shoot(1.0)))*dgamma;
  return gamma;
}


double GRBShell::beta(const double gamma) 
{
  if(gamma<1.0)
    {
      cout << "warning: gamma undefined, returning beta=0" << endl;
      return 0;
    } else 
      {
	return sqrt(1. - 1./(gamma*gamma));  
      }
}


void GRBShell::evolve(double time) 
{
  m_radius += beta(m_gamma)*cst::c*time;
}
