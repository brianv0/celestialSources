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
#include <stdexcept>

#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"

#include "facilities/Util.h"

#include "astro/SkyDir.h"

#include "flux/SpectrumFactory.h"
#include "flux/EventSource.h"

#include "FitsImage.h"
#include "Util.h"

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
   makeIntegralDistribution(m_image);

//    std::cerr << "Integral over the map: " 
//              << m_mapIntegral << std::endl;
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

void MapSource::
samplePixel(unsigned int indx, double &lon, double &lat) const {

   unsigned int i = indx % m_lon.size();
   unsigned int j = indx/m_lon.size();

// Sample uniformly in longitude
   double xi = RandFlat::shoot();
   double lon_step;
   if (i == m_lon.size()-1) {
      lon_step = m_lon.at(i) - m_lon.at(i-1);
   } else {
      lon_step = m_lon.at(i+1) - m_lon.at(i);
   }

   lon = xi*lon_step + m_lon.at(i);

// Sample as cos(lat) in latitude
   xi = RandFlat::shoot();
   double lat_step;
   if (j == m_lat.size()-1) {
      lat_step = m_lat.at(j) - m_lat.at(j-1);
   } else {
      lat_step = m_lat.at(j+1) - m_lat.at(j);
   }
   double arg = 2.*xi*cos(m_lat.at(j)*M_PI/180.)*sin(lat_step/2.*M_PI/180.)
      + sin((m_lat.at(j) - lat_step/2.)*M_PI/180.);
   lat = asin(arg)*180./M_PI;
}

double MapSource::mapValue(unsigned int i, unsigned int j) {
   unsigned int indx = j*m_lon.size() + i;
   return m_image[indx];
}

void MapSource::readFitsFile(std::string fitsFile) {
   facilities::Util::expandEnvVar(&fitsFile);

   genericSources::Util::file_ok(fitsFile);
   genericSources::FitsImage fitsImage(fitsFile);

   fitsImage.getAxisNames(m_axisTypes);

   fitsImage.getAxisVector(0, m_lon);
   fitsImage.getAxisVector(1, m_lat);

   fitsImage.getCelestialArrays(m_lonArray, m_latArray);

   fitsImage.getSolidAngles(m_solidAngles);
   fitsImage.getImageData(m_image);
}

void MapSource::
makeIntegralDistribution(const std::vector<double> & pixelValues) {
   unsigned int npix = m_solidAngles.size();
   if (pixelValues.size() < npix) {
      throw std::runtime_error("MapSource::makeIntegralDistribution:\n"
                               + std::string("pixelValues vector has fewer ")
                               + "elements than the number of image pixels");
   }
   m_integralDist.resize(npix);
   m_integralDist[0] = 0;
   double totalSolidAngle(0);
   for (unsigned int i = 1; i < npix; i++) {
      m_integralDist[i] = m_integralDist[i-1] 
         + m_solidAngles[i]*pixelValues[i];
      totalSolidAngle += m_solidAngles[i];
   }
   m_mapIntegral = m_integralDist[npix-1];
   for (unsigned int i = 1; i < npix; i++) {
      m_integralDist[i] /= m_integralDist[npix-1];
   }
//    std::cerr << "total solid angle in map: " 
//              << totalSolidAngle/M_PI << "*pi" << std::endl;
}
