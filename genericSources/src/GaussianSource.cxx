/**
 * @file GaussianSource.cxx
 * @brief A simple Spectrum subclass that exercises the flux package.
 * @author J. Chiang
 *
 * $Header$
 */

#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"

#include "facilities/Util.h"

#include "astro/SkyDir.h"

#include "flux/SpectrumFactory.h"
#include "flux/EventSource.h"

#include "genericSources/GaussianSource.h"

ISpectrumFactory &GaussianSourceFactory() {
   static SpectrumFactory<GaussianSource> myFactory;
   return myFactory;
}

GaussianSource::GaussianSource(const std::string &paramString) 
   : m_flux(1.), m_gamma(2), m_major(1), m_minor(1), m_posAngle(0),
     m_emin(30.), m_emax(1e5) {

   std::vector<std::string> params;
   facilities::Util::stringTokenize(paramString, ", ", params);

   m_flux = ::atof(params[0].c_str());
   m_gamma = ::atof(params[1].c_str());
   double ra = ::atof(params[2].c_str());
   double dec = ::atof(params[3].c_str());
   astro::SkyDir srcDir(ra, dec);
   if (params.size() > 4) m_major = ::atof(params[4].c_str());
   if (params.size() > 5) m_minor = ::atof(params[5].c_str());
   if (params.size() > 6) m_posAngle = ::atof(params[6].c_str());
   if (params.size() > 7) m_emin = ::atof(params[7].c_str());
   if (params.size() > 8) m_emax = ::atof(params[8].c_str());

// Convert to radians.
   m_major *= M_PI/180.;
   m_minor *= M_PI/180.;
   m_posAngle *= M_PI/180.;

// Rotation to align source direction with z-axis.
   m_rot = HepRotation().rotateZ(-ra*M_PI/180.).rotateY((dec - 90.)*M_PI/180.);

// Perform the rotation.
   m_rotatedSrcVec = m_rot*srcDir.dir().unit();
}

float GaussianSource::operator()(float xi) const {
   double one_m_gamma = 1. - m_gamma;
   double arg = xi*(pow(m_emax, one_m_gamma) - pow(m_emin, one_m_gamma)) 
      + pow(m_emin, one_m_gamma);
   float energy = pow(arg, 1./one_m_gamma);
   return energy;
}

double GaussianSource::flux(double time) const {
   (void)(time);
   return m_flux;
}

double GaussianSource::solidAngle() const {
   return 1;
}

double GaussianSource::interval(double time) {
   double rate = flux(time)*EventSource::totalArea();
   double xi = RandFlat::shoot();
   return -log(1. - xi)/rate;
}

double GaussianSource::energy(double time) {
   (void)(time);
   double xi = RandFlat::shoot();
   return (*this)(xi);
}

std::pair<double, double> GaussianSource::dir(double energy) {
   (void)(energy);

   double psi = RandFlat::shoot()*2.*M_PI;

   double x = m_major*cos(psi);
   double y = m_minor*sin(psi);

   double sigma = sqrt(x*x + y*y);
   double phi = atan2(y, x);

   double theta = RandGauss::shoot()*sigma;

   Hep3Vector appDir 
      = HepRotation().rotateY(theta).rotateZ(phi+m_posAngle)*m_rotatedSrcVec;

   astro::SkyDir myDir(m_rot.inverse()*appDir);
   
   return std::make_pair(myDir.l(), myDir.b());
}
