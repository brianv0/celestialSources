/**
 * @file Pulsar.cxx
 * @brief A pulsar with period, period derivative, whose light curve is
 * given by a template file.
 * @author J. Chiang
 *
 * $Header$
 */

#include <cmath>
#include <cstring>
#include <algorithm>

#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"

#include "facilities/Util.h"

#include "flux/SpectrumFactory.h"
#include "flux/EventSource.h"

#include "genericSources/Pulsar.h"

namespace {
   void readLines(std::string inputFile, 
                  std::vector<std::string> &lines) {
      std::ifstream file(inputFile.c_str());
      lines.clear();
      std::string line;
      while (std::getline(file, line, '\n')) {
         if (line != "" && line != " ") { //skip (most) blank lines
            lines.push_back(line);
         }
      }
   }
}

ISpectrumFactory &PulsarFactory() {
   static SpectrumFactory<Pulsar> myFactory;
   return myFactory;
}

Pulsar::Pulsar(const std::string &paramString)
   : PeriodicSource(2., 30., 1e5), m_meanFlux(1), m_period(0.033), 
     m_pdot(0), m_t0(0), m_phi0(0) {

   std::vector<std::string> params;
   facilities::Util::stringTokenize(paramString, ", ", params);

   m_meanFlux = std::atof(params[0].c_str());
   m_gamma = std::atof(params[1].c_str());
   m_period = std::atof(params[2].c_str());
   m_pdot = std::atof(params[3].c_str());
   m_t0 = std::atof(params[4].c_str());
   std::string templateFile = params[5];
   if (params.size() > 6) m_phi0 = std::atof(params[6].c_str());
   if (params.size() > 7) m_emin = std::atof(params[7].c_str());
   if (params.size() > 8) m_emax = std::atof(params[8].c_str());

   computeIntegralDist(templateFile);
}

double Pulsar::interval(double current_time) {
   double time = drawTime();
   double resid = residualTime(time + current_time);
   double my_period = period(time + current_time);
   double xi = resid/my_period;
   double my_interval 
      = interpolate(m_integralDist, m_phases, xi)*my_period + time - resid;
   return my_interval;
}

void Pulsar::computeIntegralDist(std::string templateFile) {
   facilities::Util::expandEnvVar(&templateFile);

// Assume the phases are bin-centered and spaced uniformly.
   std::vector<double> light_curve;
   readLightCurve(templateFile, light_curve);

   int npts = light_curve.size() + 1;
   double phase_step = 1./static_cast<double>(npts-1);
   m_phases.resize(npts);
   m_integralDist.resize(npts);

   m_phases[0] = 0.;
   m_integralDist[0] = 0.;
   for (int i = 1; i < npts; i++) {
      m_phases[i] = m_phases[i-1] + phase_step;
      m_integralDist[i] = m_integralDist[i-1] + light_curve[i-1];
   }
   for (int i = 0; i < npts; i++) {
      m_integralDist[i] /= m_integralDist[npts-1];
   }
}

void Pulsar::readLightCurve(const std::string & templateFile,
                            std::vector<double> & light_curve) {
   light_curve.clear();
   std::vector<std::string> lines;
   ::readLines(templateFile, lines);
   std::vector<std::string> tokens;
   for (std::vector<std::string>::iterator it = lines.begin();
        it != lines.end(); it++) {
      facilities::Util::stringTokenize(*it, " \t", tokens);
// Assume the second column contains the light curve.
      light_curve.push_back(std::atof(tokens[1].c_str()));
   }
}

double Pulsar::residualTime(double time) const {
   double f0 = 1./m_period;
   double f1 = -m_pdot/m_period/m_period;
   double dt = time - m_t0;
   double turns = m_phi0 + f0*dt + f1*dt*dt/2.;
   double a = f1/2.;
   double b = f0 - f1*m_t0;
   double c = m_phi0 - f0*m_t0 + f1*m_t0*m_t0/2. - std::floor(turns);
   double sgn = ((b >= 0) ? 1 : -1);
   double q = -(b + sgn*std::sqrt(b*b - 4.*a*c))/2.;
   double tn = c/q;
   return time - tn;
} 

double Pulsar::period(double time) const {
   return m_period + m_pdot*(time - m_t0);
}

double Pulsar::drawTime() const {
   double xi = RandFlat::shoot();
   double time = -std::log(xi)/m_meanFlux/EventSource::totalArea();
   return time;
}
