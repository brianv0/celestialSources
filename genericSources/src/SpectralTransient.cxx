/**
 * @file SpectralTransient.cxx
 * @brief A source that varies both spectrally and in flux with temporal
 * properties given by a template file.
 * @author J. Chiang
 *
 * $Header$
 */

#include <cmath>
#include <cstdlib>

#include <algorithm>
#include <sstream>
#include <stdexcept>

#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoisson.h"

#include "facilities/Util.h"

#include "flux/SpectrumFactory.h"
#include "flux/EventSource.h"

#include "fitsio.h"

#include "ConstParMap.h"
#include "FitsImage.h"
#include "Util.h"
#include "genericSources/EblAtten.h"
#include "genericSources/SpectralTransient.h"

std::vector<double> SpectralTransient::ModelInterval::s_energies;

ISpectrumFactory & SpectralTransientFactory() {
   static SpectrumFactory<SpectralTransient> myFactory;
   return myFactory;
}

SpectralTransient::SpectralTransient(const std::string & paramString) 
   : m_emin(20.), m_emax(2e5), m_lc(0), m_z(0), m_logParabola(0), m_tau(0) {
   std::string templateFile;
   if (paramString.find("=") == std::string::npos) {
      std::vector<std::string> params;
      facilities::Util::stringTokenize(paramString, ", ", params);
      
      m_flux = std::atof(params[0].c_str());
      m_tstart = std::atof(params[1].c_str());
      m_tstop = std::atof(params[2].c_str());
      templateFile = params[3];
      if (params.size() > 4) m_emin = std::atof(params[4].c_str());
      if (params.size() > 5) m_emax = std::atof(params[5].c_str());
      if (params.size() > 6) m_lc = std::atoi(params[6].c_str());
      if (params.size() > 7) {
         m_z = std::atof(params[7].c_str());
         m_tau = new IRB::EblAtten(IRB::Hartmann);
      }
      if (params.size() > 8) m_logParabola = std::atoi(params[8].c_str());
   } else {
      std::map<std::string, std::string> params;
      facilities::Util::keyValueTokenize(paramString, ", ", params);

      genericSources::ConstParMap parmap(params);
      m_flux = parmap.value("flux");
      m_tstart = parmap.value("tstart");
      m_tstop = parmap.value("tstop");
      templateFile = parmap["templateFile"];
      try {
         m_emin = parmap.value("emin");
      } catch (...) {
      }
      try {
         m_emax = parmap.value("emax");
      } catch (...) {
      }
      try {
         m_lc = static_cast<int>(parmap.value("lc"));
      } catch (...) {
      }
      try {
         m_z = parmap.value("z");
         m_tau = new IRB::EblAtten(IRB::Hartmann);
      } catch (...) {
      }
      try {
         m_logParabola = std::atoi(parmap["useLogParabola"].c_str());
      } catch (...) {
      }
   }

   createEvents(templateFile);
}

double SpectralTransient::interval(double time) {
   std::vector<std::pair<double, double> >::const_iterator event =
      std::upper_bound(m_events.begin(), m_events.end(), 
                       std::make_pair(time, 0), compareEventTime);
   if (event != m_events.end()) {
      m_currentEnergy = event->second;
      return event->first - time;
   } 
// There should be a better way to turn off a source than this:
   return 8.64e5;
}


void SpectralTransient::createEvents(std::string templateFile) {
   facilities::Util::expandEnvVar(&templateFile);

   genericSources::Util::file_ok(templateFile);

   if (genericSources::Util::isFitsFile(templateFile)) {
      readFitsLightCurve(templateFile);
   } else {
      readLightCurve(templateFile);
   }
   rescaleLightCurve();

   unsigned int npts(m_lightCurve.size());
   std::vector<double> tt(npts+1);
   std::vector<double> integralDist(npts+1);

   integralDist[0] = 0;

   for (unsigned int i = 0; i < npts; i++) {
      integralDist[i+1] = integralDist[i] + m_lightCurve[i].flux
         *(m_lightCurve[i].stopTime - m_lightCurve[i].startTime);
   }
   for (unsigned int i = 0; i < npts+1; i++) {
      integralDist[i] /= integralDist[npts];
   }

   double npred = m_flux*EventSource::totalArea()*(m_tstop - m_tstart);
   long nevts = RandPoisson::shoot(npred);
//    std::cerr << "SpectralTransient: number of events = "
//              << nevts << std::endl;
   m_events.clear();
   m_events.reserve(nevts);
   for (long i = 0; i < nevts; i++) {
      std::pair<double, double> my_event = drawEvent(integralDist);
      if (m_z == 0 || m_tau == 0 ||
          RandFlat::shoot() < std::exp(-(*m_tau)(my_event.second, m_z))) {
         m_events.push_back(my_event);
      }
   }
   std::stable_sort(m_events.begin(), m_events.end(), compareEventTime);
}

void SpectralTransient::
readFitsLightCurve(const std::string & templateFile) {
   std::string routineName("SpectralTransient::readFitsLightCurve");
   int hdu = genericSources::FitsImage::findHdu(templateFile, "TIMES");

   int status(0);
   fitsfile * fptr = 0;
   fits_open_file(&fptr, templateFile.c_str(), READONLY, &status);
   genericSources::FitsImage::fitsReportError(status, routineName);

   int hdutype(0);
   fits_movabs_hdu(fptr, hdu, &hdutype, &status);
   genericSources::FitsImage::fitsReportError(status, routineName);

   std::vector<double> startTimes;
   genericSources::FitsImage::readColumn(fptr, "Time", startTimes);
   unsigned int nelements(startTimes.size());

   std::vector<double> stopTimes(nelements);
   std::copy(startTimes.begin() + 1, startTimes.end(), stopTimes.begin());
   stopTimes.at(nelements-1) = 2.*startTimes.at(nelements-1) 
      - startTimes.at(nelements-2);

   hdu = genericSources::FitsImage::findHdu(templateFile, "LIGHT_CURVES");
   fits_movabs_hdu(fptr, hdu, &hdutype, &status);
   genericSources::FitsImage::fitsReportError(status, routineName);

   std::vector<double> flux, gamma1, gamma2, ebreak;
   genericSources::FitsImage::readRowVector(fptr, "Flux", m_lc,
                                            nelements, flux);
   genericSources::FitsImage::readRowVector(fptr, "Gamma1", m_lc,
                                            nelements, gamma1);
   genericSources::FitsImage::readRowVector(fptr, "Gamma2", m_lc,
                                            nelements, gamma2);
   genericSources::FitsImage::readRowVector(fptr, "Ebreak", m_lc,
                                            nelements, ebreak);
   m_lightCurve.clear();
   m_lightCurve.reserve(nelements);
   for (unsigned int i = 0; i < nelements; i++) {
      double x[] = {startTimes.at(i), stopTimes.at(i), flux.at(i),
                    gamma1.at(i), gamma2.at(i), ebreak.at(i)};
      std::vector<double> data(x, x+6);
      m_lightCurve.push_back(ModelInterval(data, m_emin, m_emax,
                                           m_logParabola));
   }
}

void SpectralTransient::
readLightCurve(const std::string & templateFile) {
   std::vector<std::string> lines;
   genericSources::Util::readLines(templateFile, lines, "#");
   std::vector<std::string>::const_iterator line;
   m_lightCurve.clear();
   m_lightCurve.reserve(lines.size());
   for (line = lines.begin(); line != lines.end(); ++line) {
      m_lightCurve.push_back(ModelInterval(*line, m_emin, m_emax,
                                           m_logParabola));
   }
}

void SpectralTransient::rescaleLightCurve() {
   double tscale = (m_tstop - m_tstart)/(m_lightCurve.back().stopTime - 
                                         m_lightCurve.front().startTime);
   std::vector<ModelInterval>::iterator interval;
   double original_startTime = m_lightCurve.front().startTime;
   for (interval = m_lightCurve.begin(); interval != m_lightCurve.end();
        ++interval) {
      interval->startTime = m_tstart + tscale*(interval->startTime -
                                               original_startTime);
      interval->stopTime = m_tstart + tscale*(interval->stopTime -
                                              original_startTime);
   }
   for (unsigned int i = 0; i < m_lightCurve.size()-1; i++) {
      if (m_lightCurve[i].stopTime > m_lightCurve[i+1].startTime) {
         std::ostringstream message;
         message << "SpectralTransient::readLightCurve:\n" 
                 << "Found stop time later than start time of "
                 << "subsequent interval:\n"
                 << "stop time: " << m_lightCurve[i-1].stopTime << "\n"
                 << "start time: " << m_lightCurve[i].startTime;
         throw std::runtime_error(message.str());
      }
   }
}   

std::pair<double, double> SpectralTransient::
drawEvent(const std::vector<double> & integralDist) const {
   double xi = RandFlat::shoot();
   std::vector<double>::const_iterator it = 
      std::upper_bound(integralDist.begin(), integralDist.end(), xi);
   int indx = it - integralDist.begin() - 1;
   double my_time = (xi - integralDist[indx])
      /(integralDist[indx+1] - integralDist[indx])
      *(m_lightCurve[indx].stopTime - m_lightCurve[indx].startTime)
      + m_lightCurve[indx].startTime;
   return std::make_pair(my_time,
                         m_lightCurve[indx].drawEnergy(m_emin, m_emax));
}

double SpectralTransient::ModelInterval::
drawEnergy(double emin, double emax) const {
   double xi = RandFlat::shoot();
   double energy;
   if (s_energies.size() > 0) {
      size_t k = std::upper_bound(m_integral.begin(), m_integral.end(), xi)
         - m_integral.begin() - 1;
      energy = (xi - m_integral.at(k))/(m_integral.at(k+1) - m_integral.at(k))
         *(s_energies.at(k+1) - s_energies.at(k)) + s_energies.at(k);
   } else {
      if (drawBelowBreak(xi)) {
         energy = genericSources::Util::drawFromPowerLaw(emin, ebreak, gamma1);
      } else {
         energy = genericSources::Util::drawFromPowerLaw(ebreak, emax, gamma2);
      }
   }
   return energy;
}

SpectralTransient::ModelInterval::
ModelInterval(const std::vector<double> & data, double emin, double emax,
              int useLogParabola) {
   startTime = data.at(0);
   stopTime = data.at(1);
   if (startTime > stopTime) {
      throw std::runtime_error("SpectralTransient::ModelInterval:\n"
                               + std::string("Poorly formed entry in FITS ")
                               + std::string("template file."));
   }
   flux = data.at(2);
   gamma1 = data.at(3);
   gamma2 = data.at(4);
   ebreak = data.at(5);
   if (useLogParabola) {
      fillCumulativeDist(emin, emax);
   } else {
      brokenPowerLawFractions(emin, emax);
   }
}

SpectralTransient::ModelInterval::ModelInterval(const std::string & line,
                                                double emin, double emax,
                                                int useLogParabola) {
   std::vector<std::string> tokens;
   facilities::Util::stringTokenize(line, ", \t", tokens);   
   if (tokens.size() != 6) {
      throw std::runtime_error("SpectralTransient::ModelInterval:\n"
                               + std::string("Poorly formed line in template ")
                               + std::string("file: ") + line);
   }
   startTime = std::atof(tokens[0].c_str());
   stopTime = std::atof(tokens[1].c_str());
   if (startTime > stopTime) {
      throw std::runtime_error("SpectralTransient::ModelInterval:\n"
                               + std::string("Poorly formed line in template ")
                               + std::string("file: ") + line);
   }
   flux = std::atof(tokens[2].c_str());
   gamma1 = std::atof(tokens[3].c_str());
   gamma2 = std::atof(tokens[4].c_str());
   ebreak = std::atof(tokens[5].c_str());
   if (useLogParabola) {
      fillCumulativeDist(emin, emax);
   } else {
      brokenPowerLawFractions(emin, emax);
   }
}

void SpectralTransient::ModelInterval::
fillCumulativeDist(double emin, double emax) {
   if (s_energies.size() == 0) {
      fillEnergies(emin, emax);
   }
   m_integral.push_back(0);
   for (size_t k = 1; k < s_energies.size(); k++) {
      const double & emin(s_energies.at(k-1));
      const double & emax(s_energies.at(k));
      m_integral.push_back(m_integral.back() + 
                           (logParabola(emin) + logParabola(emax))/2.
                           *(emax - emin));
   }
   for (size_t k = 1; k < s_energies.size(); k++) {
      m_integral.at(k) /= m_integral.back();
   }
}

double SpectralTransient::ModelInterval::
logParabola(double energy) const {
   double x = energy/ebreak;
   return std::pow(x, -(gamma1 + gamma2*std::log(x)));
}

void SpectralTransient::ModelInterval::
fillEnergies(double emin, double emax) {
   size_t nee(200);
   double estep(std::log(emax/emin)/(nee-1));
   s_energies.reserve(nee);
   for (size_t k = 0; k < nee; k++) {
      s_energies.push_back(emin*std::exp(estep*k));
   }
}

void SpectralTransient::ModelInterval::
brokenPowerLawFractions(double emin, double emax) {
   double one_m_g1 = 1. - gamma1;
   double one_m_g2 = 1. - gamma2;
   if (emin < ebreak && ebreak < emax) {
      m_lowerFraction = (1. - std::pow(emin/ebreak, one_m_g1))/one_m_g1;
      m_lowerFraction = m_lowerFraction
         /(m_lowerFraction + (std::pow(emax/ebreak, one_m_g2)-1.)/one_m_g2);
   } else if (emin > ebreak) {
      m_lowerFraction = 0;
   } else if (emax < ebreak) {
      m_lowerFraction = 1.;
   }
}
                 
