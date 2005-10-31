/**
 * @file SimpleTransient.cxx
 * @brief A flaring source defined by a single active interval.
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
#include "CLHEP/Random/RandPoisson.h"

#include "facilities/Util.h"

#include "flux/SpectrumFactory.h"
#include "flux/EventSource.h"

#include "genericSources/SimpleTransient.h"

ISpectrumFactory & SimpleTransientFactory() {
   static SpectrumFactory<SimpleTransient> myFactory;
   return myFactory;
}

SimpleTransient::SimpleTransient(const std::string & paramString) 
   : m_flux(1.), m_gamma(2), m_tstart(0), m_tstop(10), m_emin(30),
     m_emax(1e5) {

   std::vector<std::string> params;
   facilities::Util::stringTokenize(paramString, ", ", params);

   m_flux = std::atof(params[0].c_str());
   m_gamma = std::atof(params[1].c_str());
   m_tstart = std::atof(params[2].c_str());
   m_tstop = std::atof(params[3].c_str());
   if (params.size() > 4) m_emin = std::atof(params[4].c_str());
   if (params.size() > 5) m_emax = std::atof(params[5].c_str());

   createEventTimes();
}

float SimpleTransient::operator()(float xi) const {
   double one_m_gamma = 1. - m_gamma;
   double arg = xi*(pow(m_emax, one_m_gamma) - pow(m_emin, one_m_gamma)) 
      + pow(m_emin, one_m_gamma);
   float energy = pow(arg, 1./one_m_gamma);
   return energy;
}

double SimpleTransient::interval(double time) {
   std::vector<double>::const_iterator eventTime =
      std::upper_bound(m_eventTimes.begin(), m_eventTimes.end(), time);
   if (eventTime != m_eventTimes.end()) {
      return *eventTime - time;
   } 
// There should be a better way to turn off a source than this:
   return 8.64e5;
}

double SimpleTransient::energy(double time) {
   (void)(time);
   double xi = RandFlat::shoot();
   return (*this)(xi);
}

void SimpleTransient::createEventTimes() {
   double duration = m_tstop - m_tstart;
   double npred = m_flux*EventSource::totalArea()*duration;
   long nevts = RandPoisson::shoot(npred);
   std::cerr << "SimpleTransient: number of events = " 
             << nevts << std::endl;
   m_eventTimes.reserve(nevts);
   for (long i = 0; i < nevts; i++) {
      double xi = RandFlat::shoot();
      m_eventTimes.push_back(duration*xi + m_tstart);
   }
   std::stable_sort(m_eventTimes.begin(), m_eventTimes.end());
}
