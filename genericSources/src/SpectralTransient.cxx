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

#include "Util.h"
#include "genericSources/SpectralTransient.h"

ISpectrumFactory & SpectralTransientFactory() {
   static SpectrumFactory<SpectralTransient> myFactory;
   return myFactory;
}

SpectralTransient::SpectralTransient(const std::string & paramString) 
   : m_emin(20.), m_emax(2e5) {
   std::vector<std::string> params;
   facilities::Util::stringTokenize(paramString, ", ", params);

   m_flux = std::atof(params[0].c_str());
   m_tstart = std::atof(params[1].c_str());
   m_tstop = std::atof(params[2].c_str());
   std::string templateFile = params[3];
   if (params.size() > 4) m_emin = std::atof(params[4].c_str());
   if (params.size() > 5) m_emax = std::atof(params[5].c_str());

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

   readLightCurve(templateFile);

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
   std::cerr << "SpectralTransient: number of events = "
             << nevts << std::endl;
   m_events.clear();
   m_events.reserve(nevts);
   for (long i = 0; i < nevts; i++) {
      m_events.push_back(drawEvent(integralDist));
   }
   std::stable_sort(m_events.begin(), m_events.end(), compareEventTime);
}

void SpectralTransient::
readLightCurve(const std::string & templateFile) {
   std::vector<std::string> lines;
   genericSources::Util::readLines(templateFile, lines, "#");
   std::vector<std::string>::const_iterator line;
   m_lightCurve.clear();
   m_lightCurve.reserve(lines.size());
   for (line = lines.begin(); line != lines.end(); ++line) {
      m_lightCurve.push_back(ModelInterval(*line, m_emin, m_emax));
   }
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
   return std::make_pair(my_time, drawEnergy(m_lightCurve[indx]));
}

double SpectralTransient::drawEnergy(const ModelInterval & interval) const {
   double xi = RandFlat::shoot();
   double energy;
   if (interval.drawBelowBreak(xi)) {
      energy = genericSources::Util::drawFromPowerLaw(m_emin, interval.ebreak,
                                                      interval.gamma1);
   } else {
      energy = genericSources::Util::drawFromPowerLaw(interval.ebreak, m_emax,
                                                      interval.gamma2);
   }
   return energy;
}

SpectralTransient::ModelInterval::ModelInterval(const std::string & line,
                                                double emin, double emax) {
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
