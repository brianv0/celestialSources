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

namespace {
   typedef std::map<std::string, std::string> ParMap_t;
   class ConstParMap {
   public:
      ConstParMap(const ParMap_t & parmap) : m_parmap(parmap) {}
      const std::string & operator[](const std::string & name) const {
         ParMap_t::const_iterator item = m_parmap.find(name);
         if (item == m_parmap.end()) {
            throw std::runtime_error("Cannot find item named " + name);
         }
         return item->second;
      }
      size_t size() const {
         return m_parmap.size();
      }
   private:
      ParMap_t m_parmap;
   };
}

ISpectrumFactory &MapSourceFactory() {
   static SpectrumFactory<MapSource> myFactory;
   return myFactory;
}

MapSource::MapSource(const std::string & paramString) 
   : m_flux(1.), m_gamma(2), m_emin(30.), m_emax(1e5) {
   
   std::string fitsFile;
   bool createSubMap(false);
   
   if (paramString.find("=") == std::string::npos) {
// use the original params string interface
      std::vector<std::string> params;
      facilities::Util::stringTokenize(paramString, ", ", params);
      
      m_flux = std::atof(params[0].c_str());
      m_gamma = std::atof(params[1].c_str());
      fitsFile = params[2];
      if (params.size() > 3) {
         m_emin = std::atof(params[3].c_str());
      }
      if (params.size() > 4) {
         m_emax = std::atof(params[4].c_str());
      }
      if (params.size() > 5) {
         try {
            m_lonMin = std::atof(params.at(5).c_str());
            m_lonMax = std::atof(params.at(6).c_str());
            m_latMin = std::atof(params.at(7).c_str());
            m_latMax = std::atof(params.at(8).c_str());
            createSubMap = true;
         } catch (...) {
            throw std::runtime_error("Error reading sub-map bounds.\n"
                                     "There must be precisely 4 parameters "
                                     "to describe a sub-map.");
         }
      }
   } else {
// use the map version
      std::map<std::string, std::string> my_parmap;
      facilities::Util::keyValueTokenize(paramString, ", ", my_parmap);
      
      ::ConstParMap parmap(my_parmap);
      m_flux = std::atof(parmap["flux"].c_str());
      m_gamma = std::atof(parmap["gamma"].c_str());
      fitsFile = parmap["fitsFile"];
      if (parmap.size() > 3) {
         m_emin = std::atof(parmap["emin"].c_str());
      }
      if (parmap.size() > 4) {
         m_emax = std::atof(parmap["emax"].c_str());
      }
      if (parmap.size() > 5) {
         try {
            m_lonMin = std::atof(parmap["lonMin"].c_str());
            m_lonMax = std::atof(parmap["lonMax"].c_str());
            m_latMin = std::atof(parmap["latMin"].c_str());
            m_latMax = std::atof(parmap["latMax"].c_str());
            createSubMap = true;
         } catch (...) {
            throw std::runtime_error("Error reading sub-map bounds.\n"
                                     "There must be precisely 4 parameters "
                                     "to describe a sub-map.");
         }
      }
   }
   readFitsFile(fitsFile, createSubMap);
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

void MapSource::readFitsFile(std::string fitsFile, bool createSubMap) {
   facilities::Util::expandEnvVar(&fitsFile);

   genericSources::Util::file_ok(fitsFile);
   genericSources::FitsImage fitsImage(fitsFile);
   
   m_mapIntegral = fitsImage.mapIntegral();
   
   if (createSubMap) {
      getSubMapAxes(fitsImage);
      fitsImage = genericSources::FitsImage::
         sampledImage(fitsImage, m_lon, m_lat, fitsImage.coordSys());

// rescale the flux by the sub-map integral
      double new_integral = fitsImage.mapIntegral();
      m_flux *= new_integral/m_mapIntegral;
      m_mapIntegral = new_integral;
   } else {
      fitsImage.getAxisVector(0, m_lon);
      fitsImage.getAxisVector(1, m_lat);
   }

   fitsImage.getAxisNames(m_axisTypes);
   fitsImage.getCelestialArrays(m_lonArray, m_latArray);

   fitsImage.getSolidAngles(m_solidAngles);
   fitsImage.getImageData(m_image);
}

void MapSource::getSubMapAxes(const genericSources::FitsImage & fitsImage) {
   std::vector<double> axis;
   fitsImage.getAxisVector(0, axis);
   double dx = std::fabs(axis.at(1) - axis.at(0));
   m_lon.clear();
   for (size_t i = 0; i < axis.size(); i++) {
      if (m_lonMin - dx < axis.at(i) && axis.at(i) < m_lonMax + dx) {
         m_lon.push_back(axis.at(i));
      }
   }
   fitsImage.getAxisVector(1, axis);
   dx = std::fabs(axis.at(1) - axis.at(0));
   m_lat.clear();
   for (size_t j = 0; j < axis.size(); j++) {
      if (m_latMin - dx < axis.at(j) && axis.at(j) < m_latMax + dx) {
         m_lat.push_back(axis.at(j));
      }
   }
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
