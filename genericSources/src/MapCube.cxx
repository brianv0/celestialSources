/**
 * @file MapCube.cxx
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
#include "Util.h"

#include "genericSources/MapCube.h"

ISpectrumFactory &MapCubeFactory() {
   static SpectrumFactory<MapCube> myFactory;
   return myFactory;
}

MapCube::MapCube(const std::string &paramString) : MapSource() {

   std::vector<std::string> params;
   facilities::Util::stringTokenize(paramString, ", ", params);

   m_flux = std::atof(params[0].c_str());
   std::string fitsFile = params[1];

   readFitsFile(fitsFile);
   readEnergyVector(fitsFile);
   makeCumulativeSpectra();
   std::vector<double> totalCounts(m_solidAngles.size());
   for (unsigned int i = 0; i < m_solidAngles.size(); i++) {
      totalCounts[i] = m_spectra[i].back();
   }
   makeIntegralDistribution(totalCounts);

//    std::cerr << "Integral over the map: " 
//              << m_mapIntegral << std::endl;
}

float MapCube::operator()(float) const {
   double xi = RandFlat::shoot();
   std::vector<double>::const_iterator it 
      = std::upper_bound(m_integralDist.begin(), m_integralDist.end(), xi);
   unsigned int indx = it - m_integralDist.begin();

   double lon, lat;
   samplePixel(indx, lon, lat);

   if (m_axisTypes[0].find_first_of("R") == 0) { // we have Equatorial coords
      astro::SkyDir myDir(lon, lat);
      m_currentDir = std::make_pair(myDir.l(), myDir.b());
   } else { // assume Galactic coordinates
      m_currentDir = std::make_pair(lon, lat);
   }
   float energy = drawEnergy(m_spectra[indx]);
   return energy;
}

double MapCube::energy(double time) {
   (void)(time);
   double xi = RandFlat::shoot();
   return (*this)(xi);
}

double MapCube::mapValue(unsigned int i, unsigned int j, unsigned int k) {
   unsigned int indx = k*m_lon.size()*m_lat.size() + j*m_lon.size() + i;
   return m_image[indx];
}

void MapCube::readEnergyVector(const std::string & fitsFile) {

   std::string routineName("readEnergyVector");

   int hdu = genericSources::FitsImage::findHdu(fitsFile, "EBOUNDS");
   
   int status(0);
   fitsfile * fptr = 0;

   fits_open_file(&fptr, fitsFile.c_str(), READONLY, &status);
   genericSources::FitsImage::fitsReportError(status, routineName);

   int hdutype(0);
   fits_movabs_hdu(fptr, hdu, &hdutype, &status);
   genericSources::FitsImage::fitsReportError(status, routineName);

   long nrows(0);
   fits_get_num_rows(fptr, &nrows, &status);
   genericSources::FitsImage::fitsReportError(status, routineName);

   readColumn(fptr, "E_MIN", m_eMin);
   readColumn(fptr, "E_MAX", m_eMax);

   fits_close_file(fptr, &status);
   genericSources::FitsImage::fitsReportError(status, routineName);
}

void MapCube::readColumn(fitsfile * fptr, const std::string & colname,
                         std::vector<double> & coldata) const {
   std::string routineName("MapCube::readColumn");
   int status(0);
   int colnum(0);
   fits_get_colnum(fptr, CASEINSEN, const_cast<char *>(colname.c_str()),
                   &colnum, &status);
   genericSources::FitsImage::fitsReportError(status, routineName);

   long nrows(0);
   fits_get_num_rows(fptr, &nrows, &status);
   genericSources::FitsImage::fitsReportError(status, routineName);

   int anynul(0), nulval(0);
   coldata.resize(nrows);
   fits_read_col(fptr, TDOUBLE, colnum, 1, 1, nrows, &nulval, &coldata[0],
                 &anynul, &status);
   genericSources::FitsImage::fitsReportError(status, routineName);
}

void MapCube::makeCumulativeSpectra() {
   m_spectra.clear();
   m_spectra.reserve(m_solidAngles.size());

   for (unsigned int j = 0; j < m_lat.size(); j++) {
      for (unsigned int i = 0; i < m_lon.size(); i++) {
         std::vector<double> counts_spectrum(m_eMin.size());
         counts_spectrum[0] = 0;
         for (unsigned int k = 1; k < m_eMin.size(); k++) {
            counts_spectrum[k] = counts_spectrum[k-1] + mapValue(i, j, k);
         }
         m_spectra.push_back(counts_spectrum);
      }
   }
}

double MapCube::drawEnergy(const std::vector<double> & spectrum) const {
   double xi = RandFlat::shoot()*spectrum.back();
   std::vector<double>::const_iterator it 
      = std::upper_bound(spectrum.begin(), spectrum.end(), xi);
   unsigned int indx = it - spectrum.begin();
   double gamma(2.);
   double value = genericSources::Util::drawFromPowerLaw(m_eMin.at(indx),
                                                         m_eMax.at(indx),
                                                         gamma);
   return value;
}
