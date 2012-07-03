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
#include <iomanip>
#include <sstream>
#include <stdexcept>

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoisson.h"

#include "facilities/Util.h"

#include "flux/SpectrumFactory.h"
#include "flux/EventSource.h"

#include "fitsio.h"

#include "celestialSources/ConstParMap.h"

#include "FitsImage.h"
#include "Util.h"
#include "eblAtten/EblAtten.h"
#include "genericSources/SpectralTransient.h"

namespace {
   double my_RandFlat_shoot() {
      static double step(1./std::pow(2., 24));
      double xi = CLHEP::RandFlat::shoot() + CLHEP::RandFlat::shoot()*step;
      if (xi > 1) {
         xi = 1.;
      }
      return xi;
   }
}

std::vector<double> SpectralTransient::ModelInterval::s_energies;

ISpectrumFactory & SpectralTransientFactory() {
   static SpectrumFactory<SpectralTransient> myFactory;
   return myFactory;
}

SpectralTransient::SpectralTransient(const std::string & paramString) 
   : m_emin(20.), m_emax(2e5), m_lc(0), m_z(0), m_logParabola(0),
     m_specFile(false), m_tau(0), m_tauScale(1), m_particle("gamma") {
   std::string templateFile;
   std::string spectrumFile("none");
   if (paramString.find("=") == std::string::npos) {
      std::vector<std::string> params;
      facilities::Util::stringTokenize(paramString, ", ", params);
      
      m_flux = std::atof(params[0].c_str());
      m_tstart = std::atof(params[1].c_str());
      m_tstop = std::atof(params[2].c_str());
      templateFile = params[3];
      if (params.size() > 4) {
         m_emin = std::atof(params[4].c_str());
      }
      if (params.size() > 5) {
         m_emax = std::atof(params[5].c_str());
      }
      if (params.size() > 6) {
         m_lc = std::atoi(params[6].c_str());
      }
      if (params.size() > 7) {
         m_z = std::atof(params[7].c_str());
      }
      if (params.size() > 8) {
         m_logParabola = std::atoi(params[8].c_str());
      }
      if (params.size() > 9) {
         IRB::EblModel eblModel = 
            static_cast<IRB::EblModel>(std::atoi(params[9].c_str()));
         try {
            m_tau = new IRB::EblAtten(eblModel);
         } catch (std::exception & eObj) {
            std::cerr << eObj.what() << "\n"
                      << "Using default, IRB::Kneiske" << std::endl;
            m_tau = new IRB::EblAtten(IRB::Kneiske);
         }
      } else {
         m_tau = new IRB::EblAtten(IRB::Kneiske);
      }
      if (params.size() > 10) {
         spectrumFile = params[10];
      }
   } else {
      celestialSources::ConstParMap parmap(paramString);
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
         m_tau = new IRB::EblAtten(IRB::Kneiske);
      } catch (...) {
      }
      try {
         m_logParabola = std::atoi(parmap["useLogParabola"].c_str());
      } catch (...) {
      }
      try {
         IRB::EblModel eblModel = 
            static_cast<IRB::EblModel>(std::atoi(parmap["eblModel"].c_str()));
         try {
            m_tau = new IRB::EblAtten(eblModel);
         } catch (std::exception & eObj) {
            std::cerr << eObj.what() << "\n"
                      << "Using default, IRB::Kneiske" << std::endl;
            m_tau = new IRB::EblAtten(IRB::Kneiske);
         }
      } catch (...) {
         m_tau = new IRB::EblAtten(IRB::Kneiske);
      }
      try {
         m_tauScale = std::atoi(parmap["tauScale"].c_str());
      } catch (...) {
      }
      try {
         spectrumFile = parmap["spectrumFile"];
      } catch (...) {
      }
      try {
         m_particle = parmap["particle"];
      } catch (...) {
      }
   }

   readModel(templateFile);
   
   if (spectrumFile != "none") {
      readSpectrum(spectrumFile);
      m_specFile = true;
   }

   m_currentInterval = m_lightCurve.begin();
}

SpectralTransient::~SpectralTransient() {
   delete m_tau;
}

double SpectralTransient::interval(double time) {
   time -= Spectrum::startTime();
   if (m_eventCache.size() == 0) {
      fillEventCache(time);
      ++m_currentInterval;
   }
   if (!m_eventCache.empty()) {
      std::pair<double, double> thisEvent = m_eventCache.back();
      m_eventCache.pop_back();
      m_currentEnergy = thisEvent.second;
      double my_interval(thisEvent.first - time);
      if (my_interval == 0) {
         throw std::runtime_error("SpectralTransient::interval:"
                                  "zero interval computed");
      }
      return my_interval;
   }
   return 3.15e8;
}

void SpectralTransient::fillEventCache(double time) {
   if (size_t(m_currentInterval - m_lightCurve.begin()) 
       >= m_lightCurve.size()) {
      return;
   }
   for ( ; m_currentInterval != m_lightCurve.end(); ++m_currentInterval) {
      if (time < m_currentInterval->startTime) {
         time = m_currentInterval->startTime;
      }
      double dt(m_currentInterval->stopTime - time);
      double flux(m_currentInterval->flux);
      double npred(flux*EventSource::totalArea()*dt);
      int nevts;
      if (npred > 0 && (nevts = CLHEP::RandPoisson::shoot(npred)) > 0) {
         std::vector< std::pair<double, double> > my_cache;
         if (m_logParabola) {
            m_currentInterval->fillCumulativeDist(m_emin, m_emax);
         }
         for (int i = 0; i < nevts; i++) {
            double energy;
            if (m_specFile) {
               energy = drawEnergy();
            } else {
               energy = m_currentInterval->drawEnergy(m_emin, m_emax);
            }
            if (m_z == 0 || m_tau == 0 || 
                CLHEP::RandFlat::shoot() < std::exp(-m_tauScale*
                                                    (*m_tau)(energy, m_z))) {
//               double eventTime(CLHEP::RandFlat::shoot()*dt + time);
               double eventTime(::my_RandFlat_shoot()*dt + time);
               my_cache.push_back(std::make_pair(eventTime, energy));
            }
         }
         if (m_logParabola) {
            m_currentInterval->clearCumulativeDist();
         }
         if (my_cache.size() > 1) {
            std::stable_sort(my_cache.begin(), my_cache.end(),
                             compareEventTime);
         }
// Remove duplicate times
         if (my_cache.size() > 0) {
            m_eventCache.push_back(my_cache.front());
         }
         for (size_t i = 1; i < my_cache.size(); i++) {
            if (m_eventCache.back().first != my_cache.at(i).first) {
               m_eventCache.push_back(my_cache.at(i));
            }
         }
         return;
      }
   }
}

void SpectralTransient::readModel(std::string templateFile) {
   facilities::Util::expandEnvVar(&templateFile);

   genericSources::Util::file_ok(templateFile);

   if (genericSources::Util::isFitsFile(templateFile)) {
      readFitsLightCurve(templateFile);
   } else {
      readLightCurve(templateFile);
   }
   for (size_t i = 0; i < m_lightCurve.size(); i++) {
      if (m_lightCurve.at(i).flux < 0) {
//          std::ostringstream message;
//          message << "Model light curve in "
//                  << templateFile << ", for light curve number "
//                  << m_lc << " contains a negative flux value.";
//          throw std::runtime_error(message.str());
         m_lightCurve.at(i).flux = 0;
      }
   }
   rescaleLightCurve();
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

   fits_close_file(fptr, &status);
   genericSources::FitsImage::fitsReportError(status);
   
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

void SpectralTransient::readSpectrum(std::string specFile) {
   facilities::Util::expandEnvVar(&specFile);
   genericSources::Util::file_ok(specFile);

   std::vector<std::string> lines;
   genericSources::Util::readLines(specFile, lines, "#");
   std::vector<std::string>::const_iterator line;

   std::vector<double> energies;
   std::vector<double> dndes;
   for (line = lines.begin(); line != lines.end(); ++line) {
      std::vector<std::string> tokens;
      facilities::Util::stringTokenize(*line, ", \t", tokens);
      energies.push_back(std::atof(tokens[0].c_str()));
      dndes.push_back(std::atof(tokens[1].c_str()));
   }
   m_energies.clear();
   m_dnde.clear();
   m_energies.push_back(m_emin);
   m_dnde.push_back(genericSources::Util::logInterpolate(energies, dndes, 
                                                         m_emin));
   for (size_t k = 0; k < energies.size(); k++) {
      if (energies.at(k) > m_emin && energies.at(k) < m_emax) {
         m_energies.push_back(energies.at(k));
         m_dnde.push_back(dndes.at(k));
      }
   }
   m_energies.push_back(m_emax);
   m_dnde.push_back(genericSources::Util::logInterpolate(energies, dndes, 
                                                         m_emax));
   m_cumulativeDist.clear();
   m_cumulativeDist.push_back(0);
   for (size_t k = 1; k < m_energies.size(); k++) {
      m_cumulativeDist.push_back(
         genericSources::Util::powerLawIntegral(m_energies.at(k-1),
                                                m_energies.at(k),
                                                m_dnde.at(k-1),
                                                m_dnde.at(k))
         + m_cumulativeDist.back());
   }
}

double SpectralTransient::drawEnergy() const {
   double xi(CLHEP::RandFlat::shoot()*m_cumulativeDist.back());
   size_t k = std::upper_bound(m_cumulativeDist.begin(),
                               m_cumulativeDist.end(), xi) 
      - m_cumulativeDist.begin() - 1;
   double gamma = -(std::log(m_dnde.at(k+1)/m_dnde.at(k))
                    /std::log(m_energies.at(k+1)/m_energies.at(k)));
   double energy = genericSources::Util::drawFromPowerLaw(m_energies.at(k),
                                                          m_energies.at(k+1),
                                                          gamma);
   return energy;
}

void SpectralTransient::rescaleLightCurve() {
   double tscale = (m_tstop - m_tstart)/(m_lightCurve.back().stopTime - 
                                         m_lightCurve.front().startTime);
   std::vector<ModelInterval>::iterator interval;
   double original_startTime = m_lightCurve.front().startTime;
   double meanFlux(0);
   for (interval = m_lightCurve.begin(); interval != m_lightCurve.end();
        ++interval) {
      interval->startTime = m_tstart + tscale*(interval->startTime -
                                               original_startTime);
      interval->stopTime = m_tstart + tscale*(interval->stopTime -
                                              original_startTime);
      meanFlux += interval->flux*(interval->stopTime - interval->startTime);
   }

   meanFlux /= (m_tstop - m_tstart);
   double fluxScale(m_flux/meanFlux);

   for (interval = m_lightCurve.begin(); interval != m_lightCurve.end();
        ++interval) {
      interval->flux *= fluxScale;
   }

   for (unsigned int i = 0; i < m_lightCurve.size()-1; i++) {
      if (m_lightCurve[i].stopTime > m_lightCurve[i+1].startTime) {
         std::ostringstream message;
         message << "SpectralTransient::rescaleLightCurve:\n" 
                 << "Found stop time later than start time of "
                 << "subsequent interval:\n"
                 << "stop time: " << m_lightCurve[i-1].stopTime << "\n"
                 << "start time: " << m_lightCurve[i].startTime;
         throw std::runtime_error(message.str());
      }
   }
}   

double SpectralTransient::ModelInterval::
drawEnergy(double emin, double emax) const {
   double xi = CLHEP::RandFlat::shoot();
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
// Call this inside of SpectralTransient::fillEventCache just before
// filling the cache
//      fillCumulativeDist(emin, emax);
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
// Call this inside of SpectralTransient::fillEventCache just before
// filling the cache
//       fillCumulativeDist(emin, emax);
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
   double lowerIntegral;
   double upperIntegral;
   if (emin < ebreak && ebreak < emax) {
      if (gamma1 == 1) {
         lowerIntegral = std::log(ebreak/emin);
      } else {
         double one_m_g1(1. - gamma1);
         lowerIntegral = (1. - std::pow(emin/ebreak, one_m_g1))/one_m_g1;
      }
      if (gamma2 == 1) {
         upperIntegral = std::log(emax/ebreak);
      } else {
         double one_m_g2(1. - gamma2);
         upperIntegral = (std::pow(emax/ebreak, one_m_g2) - 1.)/one_m_g2;
      }
      m_lowerFraction = lowerIntegral/(lowerIntegral + upperIntegral);
   } else if (emin > ebreak) {
      m_lowerFraction = 0;
   } else if (emax < ebreak) {
      m_lowerFraction = 1.;
   }
}
                 
