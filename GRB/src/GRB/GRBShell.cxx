//
//    GRBShell: Class that describes the a Shell
//    Authors: Nicola Omodei & Johann Cohen Tanugi 
//

#include "GRBShell.h"
#include <iostream>
#include <math.h>
#include <stdlib.h>

using namespace std;

GRBShell::GRBShell(double gamma, double mass,
		   double thickness, double radius, double distance,string type) 
  : m_gamma(gamma), m_mass(mass),m_thickness(thickness),m_radius(radius),m_distance(distance)
{ 
  if(type =="jet")
    {
      std::cout<<" Sell type = jet "<<std::endl;
      m_volume = cst::pi*pow(m_radius,2.)*m_thickness;
    }
  else
    {
      std::cout<<" Sell type = isotropic "<<std::endl;
      m_volume = 4.*cst::pi*pow(m_distance,2.)*m_thickness;
    }
}



double GRBShell::beta(const double gamma) 
{
  if(gamma<1.0)
    {
      return 0;
    }
  else 
    {
      return sqrt(1. - 1./(gamma*gamma));  
    }
}


void GRBShell::evolve(double dt) 
{
  //Interaction with the Inter Stellar Medium: 
  if(m_distance>1.0e+19)
    {
      m_gamma=cst::viscosity+m_gamma*(1.0-cst::viscosity);
    }
  if(m_gamma<1.0) m_gamma=1.0;
  m_distance += beta(m_gamma)*cst::c*dt;
  // Expanding sells...
  //  m_thickness=m_distance/pow(m_gamma,2)>m_thickness?m_distance/pow(m_gamma,2):m_thickness;
  
}


GRBShell GRBShell::operator+(GRBShell Sh) 
{ 
  m_distance = Sh.getDistance();
  //Kinematics:  
  double m1 = m_mass;
  double m2 = Sh.getMass();
  double g1 = m_gamma;
  double g2 = Sh.getGamma();
  double dr1= m_thickness;
  double dr2= Sh.getThickness();
  
  //Number and density of particles in the shocked material:
  //std::cout<<m_gamma<<std::endl;
  
  // Gamma of the resulting shell: 
  m_gamma=sqrt((m1*g1+m2*g2)/(m1/g1+m2/g2));
  m_gamma=((m_gamma<=1.) ? 1.+1.0e-6 : m_gamma);
  // See Piran 1999
  // Internal Energy:
  //Updating Shell Info. after shock:
  setMass(m1+m2);
  setGamma(m_gamma);
  //std::cout<<m_gamma<<std::endl;
  setThickness((dr1+dr2)/2.);
  setEint(m1*cst::c2*(g1-m_gamma)+m2*cst::c2*(g2-m_gamma));
  return (*this);
}

