/**
 * @file TransientTempate.h
 * @brief A flaring source whose light curve shape is given by a
 * template file.
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

#include "genericSources/Pulsar.h"
#include "genericSources/TransientTemplate.h"

ISpectrumFactory & TransientTemplateFactory() {
   static SpectrumFactory<TransientTemplate> myFactory;
   return myFactory;
}

TransientTemplate::TransientTemplate(const std::string & paramString) 
   : SimpleTransient() {

   std::vector<std::string> params;
   facilities::Util::stringTokenize(paramString, ", ", params);

   m_flux = std::atof(params[0].c_str());
   m_gamma = std::atof(params[1].c_str());
   m_tstart = std::atof(params[2].c_str());
   m_tstop = std::atof(params[3].c_str());
   std::string templateFile = params[4];
   if (params.size() > 5) m_emin = std::atof(params[5].c_str());
   if (params.size() > 6) m_emax = std::atof(params[6].c_str());

   createEventTimes(templateFile);
}

void TransientTemplate::createEventTimes(std::string templateFile) {
   facilities::Util::expandEnvVar(&templateFile);

   std::vector<double> light_curve;
   Pulsar::readLightCurve(templateFile, light_curve);

   unsigned int npts = light_curve.size() + 1;
   double duration = m_tstop - m_tstart;
   double tstep = duration/static_cast<double>(npts - 1);
   std::vector<double> tt(npts);
   std::vector<double> integralDist(npts);

   tt[0] = m_tstart;
   integralDist[0] = 0;
   for (unsigned int i = 1; i < npts; i++) {
      tt[i] = tt[i-1] + tstep;
      integralDist[i] = integralDist[i-1] + light_curve[i-1];
   }
   for (unsigned int i = 0; i < npts; i++) {
      integralDist[i] /= integralDist[npts-1];
   }

   double npred = m_flux*EventSource::totalArea()*duration;
   long nevts = RandPoisson::shoot(npred);
   std::cerr << "TemplateTransient: number of events = "
             << nevts << std::endl;
   m_eventTimes.reserve(nevts);
   for (long i = 0; i < nevts; i++) {
      m_eventTimes.push_back(drawTime(tt, integralDist));
   }
   std::stable_sort(m_eventTimes.begin(), m_eventTimes.end());
}

double TransientTemplate::drawTime(const std::vector<double> & tt,
                                   const std::vector<double> & integralDist) {
   double xi = RandFlat::shoot();
   std::vector<double>::const_iterator it = 
      std::upper_bound(integralDist.begin(), integralDist.end(), xi);
   int indx = it - integralDist.begin() - 1;
/// @bug Using iterators causes crash if there are zero value entries
/// in the flux light curve.
//    double my_time = (xi - *it)/(*(it+1) - *it)*(tt[indx+1] - tt[indx])
//       + tt[indx];
   double my_time = (xi - integralDist[indx])
      /(integralDist[indx+1] - integralDist[indx])*(tt[indx+1] - tt[indx])
      + tt[indx];
   return my_time;
}
