/**
 * @file microQuasar.cxx
 * @brief A phenomenological model of the microQuasar based on EGRET measurements
 * @author D. Petry
 *
 * $Header$
 */

#include <iostream>

#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <stdexcept>

#include "facilities/Util.h"

#include "flux/SpectrumFactory.h"
#include "flux/EventSource.h"
#include "CLHEP/Random/RandFlat.h"

#include "microQuasar/microQuasar.h"



ISpectrumFactory &microQuasarFactory() { // a.k.a. microQuasarFactory, see http://www.bbc.co.uk/dna/h2g2/A105265
   static SpectrumFactory<microQuasar> myFactory;
   return myFactory;
}

microQuasar::microQuasar(const std::string &paramString) {

   std::vector<std::string> params;
   facilities::Util::stringTokenize(paramString, ", ", params);

   m_ftot = 9.3e-2;
   m_alt = 565e3;
   m_eMin = 20.;;
   m_eMax = 200000.;

   m_orbitalPeriod = 3.91*86400;
   m_orbitalModulation = 1.0;
   m_phi0 = 0.;

   std::cerr << "microQuasar created. Total flux = " 
	     << m_ftot << " cm^-2 s^-1 " << " between " 
	     << m_eMin << " MeV and "
	     << m_eMax << " MeV." << std::endl;

}


float microQuasar::operator()(float xi) const {
   double one_m_gamma = 1. - m_currentSpectralIndex;
   double arg = xi*(pow(m_eMax, one_m_gamma) - pow(m_eMin, one_m_gamma)) 
      + pow(m_eMin, one_m_gamma);
   float energy = pow(arg, 1./one_m_gamma);
   return energy;
}

double microQuasar::energy(double time) {
	int region = m_orbitalRegion.findRegion(time,m_orbitalPeriod);
	m_currentSpectralIndex = m_orbitalRegion.getSpectralIndex(region);
	double xi = CLHEP::RandFlat::shoot();
   return (*this)(xi);
}


void microQuasar::modulation(const double x, double& funcValue, double& derivValue) {
// see http://d0.phys.washington.edu/~burnett/glast/generate_periodic/oscilations.htm for details
	double scale = m_ftot*EventSource::totalArea();
	double now = fmod((float)m_currentTime,(float)m_orbitalPeriod)/m_orbitalPeriod*2.*M_PI;
	float z = -log(1.-m_randPhase);
	float m_nTurns = floor(z/(2.*M_PI*scale));
	float zp = z - 2.*M_PI*m_nTurns*scale;
	funcValue =  scale*(x - m_orbitalModulation*(sin(x+now)-sin(now))) - zp;
	derivValue = scale*(1.- m_orbitalModulation*cos(x+now));
	return;
}

double microQuasar::interval(double current_time) {
	m_currentTime = current_time - Spectrum::startTime();
	m_randPhase = CLHEP::RandFlat::shoot();
	double deltaT = m_orbitalPeriod/(2.*M_PI)*(rtsafe(0.,2.*M_PI,1.e-2)+2.*M_PI*m_nTurns);
//	std::cout << " current t " << m_currentTime << "deltaT " << deltaT << std::endl;
	return deltaT;
//	return 20.;
}


int microQuasar::OrbitalRegion::findRegion(double time, float period) {
	float timef = time;
	float phase = fmod(timef,period)/period;
	return (phase > m_minOrbitalPhase && phase < m_maxOrbitalPhase) ? 1 : 2;
}
double microQuasar::rtsafe(const double x1, const double x2,	const double xacc)
{
	const int MAXIT=100;
	int j;
	double df,dx,dxold,f,fh,fl,temp,xh,xl,rts;

	modulation(x1,fl,df);
	modulation(x2,fh,df);
	if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
		std::cerr << "Root must be bracketed in rtsafe" << std::endl;
	if (fl == 0.0) return x1;
	if (fh == 0.0) return x2;
	if (fl < 0.0) {
		xl=x1;
		xh=x2;
	} else {
		xh=x1;
		xl=x2;
	}
	rts=0.5*(x1+x2);
	dxold=fabs(x2-x1);
	dx=dxold;
	modulation(rts,f,df);
	for (j=0;j<MAXIT;j++) {
		if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0)
			|| (fabs(2.0*f) > fabs(dxold*df))) {
			dxold=dx;
			dx=0.5*(xh-xl);
			rts=xl+dx;
			if (xl == rts) return rts;
		} else {
			dxold=dx;
			dx=f/df;
			temp=rts;
			rts -= dx;
			if (temp == rts) return rts;
		}
		if (fabs(dx) < xacc) return rts;
		modulation(rts,f,df);
		if (f < 0.0)
			xl=rts;
		else
			xh=rts;
	}
	std::cerr << "Maximum number of iterations exceeded in rtsafe" << std::endl;
	return 0.0;
}

