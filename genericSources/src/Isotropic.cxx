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

#include "facilities/Util.h"

#include "flux/SpectrumFactory.h"
#include "flux/EventSource.h"

#include "genericSources/Isotropic.h"

ISpectrumFactory &IsotropicFactory() {
   static SpectrumFactory<Isotropic> myFactory;
   return myFactory;
}

Isotropic::Isotropic(const std::string & paramString) 
   : m_flux(1.), m_gamma(2), m_emin(30.), m_emax(1e5) {
   std::vector<std::string> params;
   facilities::Util::stringTokenize(paramString, ", ", params);

   m_flux = std::atof(params[0].c_str());
   m_gamma = std::atof(params[1].c_str());
   if (params.size() > 2) m_emin = std::atof(params[2].c_str());
   if (params.size() > 3) m_emax = std::atof(params[3].c_str());
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
   double xi = RandFlat::shoot();
   return -log(1. - xi)/rate;
}

double Isotropic::energy(double time) {
   (void)(time);
   double xi = RandFlat::shoot();
   return (*this)(xi);
}

std::pair<double, double> Isotropic::dir(double energy) {
   (void)(energy);

   double xi = RandFlat::shoot();
   double lon = 360.*xi - 180.;

   xi = RandFlat::shoot();
   double mu = 2.*xi - 1.;
   double theta;
   if (fabs(mu) <= 1.) {
      theta = acos(mu);
   } else {
      if (mu > 0) {
         theta = 0;
      } else {
         theta = M_PI;
      }
   }
   double lat = (M_PI/2. - theta)*180./M_PI;

   return std::make_pair(lon, lat);
}
