/**
 * @file MapSource.cxx
 * @brief A simple Spectrum subclass that exercises the flux package.
 * @author J. Chiang
 *
 * $Header$
 */

#include <cmath>
#include <cstdlib>

#include <algorithm>

#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"

#include "facilities/Util.h"

#include "astro/SkyDir.h"

#include "flux/SpectrumFactory.h"
#include "flux/EventSource.h"

#include "FitsImage.h"
#include "genericSources/MapSource.h"

ISpectrumFactory &MapSourceFactory() {
   static SpectrumFactory<MapSource> myFactory;
   return myFactory;
}

MapSource::MapSource(const std::string &paramString) 
   : m_flux(1.), m_gamma(2), m_emin(30.), m_emax(1e5) {

   std::vector<std::string> params;
   facilities::Util::stringTokenize(paramString, ", ", params);

   m_flux = std::atof(params[0].c_str());
   m_gamma = std::atof(params[1].c_str());
   std::string fitsFile = params[2];
   if (params.size() > 3) m_emin = std::atof(params[3].c_str());
   if (params.size() > 4) m_emax = std::atof(params[4].c_str());

   readFitsFile(fitsFile);

   std::cerr << "Integral over the map: " 
             << m_mapIntegral << std::endl;
}

float MapSource::operator()(float xi) const {
   double one_m_gamma = 1. - m_gamma;
   double arg = xi*(pow(m_emax, one_m_gamma) - pow(m_emin, one_m_gamma)) 
      + pow(m_emin, one_m_gamma);
   float energy = pow(arg, 1./one_m_gamma);
   return energy;
}

double MapSource::flux(double time) const {
   (void)(time);
   return m_flux;
}

double MapSource::solidAngle() const {
   return 1;
}

double MapSource::interval(double time) {
   double rate = flux(time)*EventSource::totalArea();
   double xi = RandFlat::shoot();
   return -log(1. - xi)/rate;
}

double MapSource::energy(double time) {
   (void)(time);
   double xi = RandFlat::shoot();
   return (*this)(xi);
}

std::pair<double, double> MapSource::dir(double energy) {
   (void)(energy);

   double xi = RandFlat::shoot();
   std::vector<double>::const_iterator it 
      = std::upper_bound(m_integralDist.begin(), m_integralDist.end(), xi);
   unsigned int indx = it - m_integralDist.begin();

   double lon, lat;
   samplePixel(indx, lon, lat);

   if (m_axisTypes[0].find_first_of("R") == 0) {
// We have Equatorial coordinates.
      astro::SkyDir myDir(lon, lat);
      return std::make_pair(myDir.l(), myDir.b());
   }
// Assume Galactic coordinates by default.
   return std::make_pair(lon, lat);
}

void MapSource::samplePixel(unsigned int indx, double &lon, double &lat) {

   unsigned int i = indx % m_lon.size();
   unsigned int j = indx/m_lon.size();

// Sample uniformly in longitude
   double xi = RandFlat::shoot();
   double lon_step = m_lon[i+1] - m_lon[i];
   lon = xi*lon_step + m_lon[i];

// Sample as cos(lat) in latitude
   xi = RandFlat::shoot();
   double lat_step = m_lat[j+1] - m_lat[j];
   double arg = 2.*xi*cos(m_lat[j]*M_PI/180.)*sin(lat_step/2.*M_PI/180.)
      + sin((m_lat[j] - lat_step/2.)*M_PI/180.);
   lat = asin(arg)*180./M_PI;
}

double MapSource::mapValue(unsigned int i, unsigned int j) {
   unsigned int indx = j*m_lon.size() + i;
   return m_image[indx];
}

void MapSource::readFitsFile(std::string fitsFile) {
   facilities::Util::expandEnvVar(&fitsFile);

   genericSources::FitsImage fitsImage(fitsFile);

   fitsImage.fetchAxisNames(m_axisTypes);

   fitsImage.fetchAxisVector(0, m_lon);
   fitsImage.fetchAxisVector(1, m_lat);

   fitsImage.fetchCelestialArrays(m_lonArray, m_latArray);

   fitsImage.fetchSolidAngles(m_solidAngles);
   fitsImage.fetchImageData(m_image);

   int npix = m_solidAngles.size();
   m_integralDist.resize(npix);
   m_integralDist[0] = 0;
   for (int i = 1; i < npix; i++) {
      m_integralDist[i] = m_integralDist[i-1] + m_solidAngles[i]*m_image[i];
   }
   m_mapIntegral = m_integralDist[npix-1];
   for (int i = 1; i < npix; i++) {
      m_integralDist[i] /= m_integralDist[npix-1];
   }
}
