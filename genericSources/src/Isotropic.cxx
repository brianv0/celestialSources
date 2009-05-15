/**
 * @file Isotropic.cxx
 * @brief Isotropic diffuse emission
 * @author J. Chiang
 *
 * $Header$
 */

#include <cmath>
#include <cstdlib>

#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "facilities/Util.h"

#include "astro/SkyDir.h"

#include "flux/SpectrumFactory.h"
#include "flux/EventSource.h"

#include "genericSources/Isotropic.h"

#include "ConstParMap.h"

ISpectrumFactory &IsotropicFactory() {
   static SpectrumFactory<Isotropic> myFactory;
   return myFactory;
}

Isotropic::Isotropic(const std::string & paramString) 
   : m_flux(1.), m_gamma(2), m_emin(20.), m_emax(5e5), m_ra(0), m_dec(0),
     m_cos_thetamax(-1) {

   if (paramString.find("=") == std::string::npos) {
      std::vector<std::string> params;
      facilities::Util::stringTokenize(paramString, ", ", params);
      
      m_flux = std::atof(params[0].c_str());
      m_gamma = std::atof(params[1].c_str());
      m_emin = std::atof(params[2].c_str());
      m_emax = std::atof(params[3].c_str());
      if (params.size() > 4) {
         m_ra = std::atof(params[4].c_str());
         m_dec = std::atof(params[5].c_str());
         m_cos_thetamax = std::cos(std::atof(params[6].c_str())*M_PI/180.);
      }
   } else {
      genericSources::ConstParMap parmap(paramString);
      
      m_flux = parmap.value("flux");
      m_gamma = parmap.value("gamma");
      m_emin = parmap.value("emin");
      m_emax = parmap.value("emax");
      try {
         m_ra = parmap.value("ra");
         m_dec = parmap.value("dec");
         m_cos_thetamax = std::cos(parmap.value("radius")*M_PI/180.);
      } catch (...) {
      }
   }
}

float Isotropic::operator()(float xi) const {
   double one_m_gamma = 1. - m_gamma;
   double arg = xi*(pow(m_emax, one_m_gamma) - pow(m_emin, one_m_gamma)) 
      + pow(m_emin, one_m_gamma);
   float energy = pow(arg, 1./one_m_gamma);
   return energy;
}

double Isotropic::flux(double time) const {
   (void)(time);
   return m_flux;
}

double Isotropic::solidAngle() const {
   return 1.;
}

double Isotropic::interval(double time) {
   double rate = flux(time)*EventSource::totalArea();
   double xi = CLHEP::RandFlat::shoot();
   return -log(1. - xi)/rate;
}

double Isotropic::energy(double time) {
   (void)(time);
   double xi = CLHEP::RandFlat::shoot();
   return (*this)(xi);
}

std::pair<double, double> Isotropic::dir(double energy) {
   (void)(energy);

   double xi = CLHEP::RandFlat::shoot();
   double phi = 2.*M_PI*xi;

   xi = CLHEP::RandFlat::shoot();
   double costh = 1. - xi*(1. - m_cos_thetamax);
   double sinth = std::sqrt(1. - costh*costh);

   CLHEP::Hep3Vector direction(std::cos(phi)*sinth, 
                               std::sin(phi)*sinth, 
                               costh);
                               
   astro::SkyDir dir(direction.rotateY((90 - m_dec)*M_PI/180.).rotateZ(m_ra*M_PI/180));

   return std::make_pair(dir.l(), dir.b());
}
