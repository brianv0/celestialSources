/**
 * @file PeriodicSource.cxx
 * @brief A simple Spectrum subclass that exercises the flux package.
 * @author J. Chiang
 *
 * $Header$
 */

#include <cmath>
#include <algorithm>

#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"

#include "facilities/Util.h"

#include "flux/SpectrumFactory.h"
#include "flux/EventSource.h"

#include "genericSources/PeriodicSource.h"

ISpectrumFactory &PeriodicSourceFactory() {
   static SpectrumFactory<PeriodicSource> myFactory;
   return myFactory;
}

PeriodicSource::PeriodicSource(const std::string &paramString) 
   : m_gamma(2), m_emin(30.), m_emax(1e5), m_flux(1.), m_period(1e3), 
     m_amplitude(0.5), m_phi0(0) {

   std::vector<std::string> params;
   facilities::Util::stringTokenize(paramString, ", ", params);

   m_flux = ::atof(params[0].c_str());
   m_gamma = ::atof(params[1].c_str());
   m_period = ::atof(params[2].c_str());
   if (params.size() > 3) m_amplitude = ::atof(params[3].c_str());
   if (params.size() > 4) m_phi0 = ::atof(params[4].c_str());
   if (params.size() > 5) m_emin = ::atof(params[5].c_str());
   if (params.size() > 6) m_emax = ::atof(params[6].c_str());

// Convert to radians from unit interval.
   m_phi0 *= 2.*M_PI;

   computeIntegralDistribution();
}

float PeriodicSource::operator()(float xi) const {
   double one_m_gamma = 1. - m_gamma;
   double arg = xi*(pow(m_emax, one_m_gamma) - pow(m_emin, one_m_gamma)) 
      + pow(m_emin, one_m_gamma);
   float energy = pow(arg, 1./one_m_gamma);
   return energy;
}

double PeriodicSource::energy(double time) {
   (void)(time);
   double xi = RandFlat::shoot();
   return (*this)(xi);
}

double PeriodicSource::interval(double time) {
   double phase = fmod(time, m_period);
   std::vector<double>::const_iterator it 
      = std::upper_bound(m_arrTimes.begin(), m_arrTimes.end(),
                         phase) - 1;
   unsigned int imin = it - m_arrTimes.begin();
   unsigned int npts = m_arrTimes.size();
   std::vector<double> newDist(npts);
   std::transform(m_integralDist.begin() + imin,
                  m_integralDist.begin() + imin +  npts,
                  newDist.begin(), 
                  std::bind2nd(std::minus<double>(), m_integralDist[imin]));
   double xi = -log(RandFlat::shoot());
   unsigned int turns = static_cast<unsigned int>(xi/newDist[npts-1]);
   double resid = fmod(xi, newDist[npts-1]);
   double my_interval = interpolate(newDist, m_arrTimes, resid) 
      + turns*m_period;
   return my_interval;
}

double PeriodicSource::interpolate(const std::vector<double> &x,
                                   const std::vector<double> &y,
                                   double xx) {
   std::vector<double>::const_iterator it 
      = std::upper_bound(x.begin(), x.end(), xx) - 1;
   unsigned int indx = it - x.begin();
   double yy;
/// @bug Using iterators causes crash if there are zero value entries
/// in the flux light curve.
//    if (*(it+1) != *it) {
//       yy = (xx - *it)/(*(it+1) - *it)*(y[indx+1] - y[indx]) + y[indx];
//    } else {
//       yy = (y[indx+1] + y[indx])/2.;
//    }
   if (x[indx+1] != x[indx]) {
      yy = (xx - x[indx])/(x[indx+1] - x[indx])
         *(y[indx+1] - y[indx]) + y[indx];
   } else {
      yy = (y[indx+1] + y[indx])/2.;
   }
   return yy;
}

void PeriodicSource::computeIntegralDistribution() {
   static unsigned int npts = 1000;
   makeGrid(npts, 0., m_period, m_arrTimes);

   double angularFreq = 2.*M_PI/m_period;

   m_integralDist.clear();
   m_integralDist.reserve(2*npts);
   m_integralDist.push_back(0);
   std::vector<double>::const_iterator it = m_arrTimes.begin();
   for ( ; it != m_arrTimes.end(); it++) {
      m_integralDist.push_back(m_flux*(*it + m_amplitude/angularFreq
                                       *(cos(angularFreq*(*it) + m_phi0)
                                         - cos(m_phi0))));
   }

// Extend to cover two periods.
   for (unsigned int i = 0; i < npts; i++) {
      m_integralDist.push_back(m_integralDist[i] + m_integralDist[npts-1]);
   }
}

void PeriodicSource::makeGrid(unsigned int n, double xmin, double xmax,
                              std::vector<double> &x, bool makeLog) {
   double xstep;
   if (makeLog) {
      xstep = log(xmax/xmin)/(n-1);
   } else {
      xstep = (xmax - xmin)/(n-1);
   }

   x.clear();
   x.reserve(n);
   for (unsigned int i = 0; i < n; i++) {
      if (makeLog) {
         x.push_back(xmin*exp(xstep*i));
      } else {
         x.push_back(xstep*i + xmin);
      }
   }
}
