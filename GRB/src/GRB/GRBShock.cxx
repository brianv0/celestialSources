//
//    GRBShock: Class that describes the a Shock between 2 Shells
//    Authors: Nicola Omodei & Johann Cohen Tanugi 
//

#include <iostream>
#include <math.h>
#include "GRBShock.h"

using namespace std;

double Bfield(double Ub)
{
  // Ub in erg, Bfield in Gauss
  return 1.0e+4*sqrt(200.0*M_PI*cst::mu0*Ub);
}

GRBShock::GRBShock(GRBShell Sh1) 
{
  m_mass      = Sh1.getMass();
   
  m_thickness = Sh1.getThickness();
  m_volume    = Sh1.getVolCom(); // [cm^3]
  m_gf        = Sh1.getGamma();
  m_Eint      = Sh1.getEint();
  
  m_partnumber  = m_mass*cst::c2/(cst::mpc2)*cst::erg2MeV;
  m_partdensity = m_partnumber/m_volume;
  m_Ue          = cst::alphae*m_Eint/m_volume;
  m_Ub          = cst::alphab*m_Eint/m_volume;
  m_Beq         = Bfield(m_Ub);
  // Acceleration of the electrons:
  double c1 = 1.e+5; 
  m_gemin = (cst::p-2.)/(cst::p-1.)*(m_Ue*cst::erg2MeV)/(m_partdensity*cst::mec2)
    *(1.-pow(c1,1.-cst::p))/(1.-pow(c1,2.-cst::p));
  m_gemin=(m_gemin <= 1.) ? 1.+1.0e-7: m_gemin;
  m_gemax = c1*m_gemin;
}

//////////////////////////////////////////////////
// Syncrotron Radiation:
//////////////////////////////////////////////////

double GRBShock::duration()
{
  double tau_rise    = m_thickness/(cst::c)/m_gf;
  double tau_decay   = 1.e+7*m_volume/m_Eint/m_gf;
  //  std::cout<<" Rise  Time "<<tau_rise<<std::endl;
  //  std::cout<<" Decay Time "<<tau_decay<<std::endl;
  return (1./2.7*sqrt(tau_decay)+tau_rise);
}
/*
  std::vector<double> GRBShock::FluxAtT(double time)
  {
  std::vector<double> spectrum(cst::enstep,0.);
  if (time<=m_tobs) return spectrum;
  GRBSynchrotron m_syn(time-m_tobs,m_gf,m_Beq,
  m_partnumber,m_gemin,m_gemax,
  m_volume,m_thickness);
  //std::cout<<m_syn->getSpectrum().size()<<std::endl;
  return m_syn.getSpectrum();
  }
*/

void GRBShock::Write()
{
  std::cout<< "--------------------Shock's parameters --------------" << std::endl;
  std::cout<< " T obs                      = "<<m_tobs<<std::endl;
  std::cout<< " Mass [g]                   = "<<m_mass<<std::endl;
  std::cout<< " Thickness [cm]             = "<< m_thickness<<std::endl;
  std::cout<< " Volume [cm^3]              = "<< m_volume<<std::endl;
  std::cout<< " Gamma final                = "<<m_gf<<std::endl;
  std::cout<< " Internal Energy [erg]      = "<<m_Eint<<std::endl;
  std::cout<< " Electron density [1/cm^3]  = "<<m_partdensity<<std::endl;
  std::cout<< " Number of Electron         = "<<m_partnumber<<std::endl;
  std::cout<< " Electron minimum gamma     = "<< m_gemin<<std::endl;
  std::cout<< " Electron Maximum gamma     = "<< m_gemax<<std::endl;
  std::cout<< " Electron Energy [erg/cm^3] = "<< m_Ue<<std::endl;
  std::cout<< " Magnetic Field [Gauss]     = "<< m_Beq<<std::endl;
  std::cout<< " Magnetic Energy [erg/cm^3] = "<< m_Ub<<std::endl;
  std::cout<< " Rise  time                 = "<<m_thickness/(cst::c)/m_gf<<std::endl;
  std::cout<< " Decay time                 = "<<1.e+7*m_volume/m_Eint/m_gf<<std::endl;
  std::cout<< " Duration                   = "<<duration()<<std::endl;
}

