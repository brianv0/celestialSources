/**
 * @file IsotropicFileSpectrum.cxx
 * @brief IsotropicFileSpectrum diffuse emission
 * @author J. Chiang
 *
 * $Header$
 */

#include <cmath>
#include <cstdlib>

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "facilities/Util.h"

#include "astro/SkyDir.h"

#include "flux/EventSource.h"
#include "flux/SpectrumFactory.h"

#include "genericSources/FileSpectrum.h"
#include "genericSources/IsotropicFileSpectrum.h"

#include "celestialSources/ConstParMap.h"

ISpectrumFactory & IsotropicFileSpectrumFactory() {
   static SpectrumFactory<IsotropicFileSpectrum> myFactory;
   return myFactory;
}

IsotropicFileSpectrum::IsotropicFileSpectrum(const std::string & paramString) 
   : m_fileSpectrum(new FileSpectrum(paramString)),
     m_ra(0), m_dec(0), m_cos_thetamax(-1) {
   celestialSources::ConstParMap parmap(paramString);
   try {
      m_ra = parmap.value("ra");
      m_dec = parmap.value("dec");
      m_cos_thetamax = std::cos(parmap.value("radius")*M_PI/180.);
   } catch (...) {
   }
}

float IsotropicFileSpectrum::operator()(float xi) const {
   return m_fileSpectrum->operator()(xi);
}

double IsotropicFileSpectrum::flux(double time) const {
   (void)(time);
   return m_fileSpectrum->flux();
}

double IsotropicFileSpectrum::solidAngle() const {
   return 1.;
}

double IsotropicFileSpectrum::interval(double time) {
   double rate = flux(time)*EventSource::totalArea();
   double xi = CLHEP::RandFlat::shoot();
   return -log(1. - xi)/rate;
}

double IsotropicFileSpectrum::energy(double time) {
   (void)(time);
   double xi = CLHEP::RandFlat::shoot();
   return (*this)(xi);
}

std::pair<double, double> IsotropicFileSpectrum::dir(double energy) {
   (void)(energy);

   double xi = CLHEP::RandFlat::shoot();
   double phi = 2.*M_PI*xi;

   xi = CLHEP::RandFlat::shoot(); 
   double costh = 1. - xi*(1. - m_cos_thetamax);
   double sinth = std::sqrt(1. - costh*costh);

   CLHEP::Hep3Vector direction(std::cos(phi)*sinth, 
                               std::sin(phi)*sinth, 
                               costh);
                               
   astro::SkyDir dir(direction.rotateY((90 - m_dec)*M_PI/180.)
                     .rotateZ(m_ra*M_PI/180));

   return std::make_pair(dir.l(), dir.b());
}
