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

microQuasar::microQuasar(const std::string &paramString) 
: m_currentTime(0)
, m_nTurns(0)
, m_currentSpectralIndex(2)
, m_randPhase(0)
{

	std::vector<std::string> params;
	facilities::Util::stringTokenize(paramString, ", ", params);

	float daySecs = 86400.;

	m_ftot = ::atof(params[0].c_str());
	m_eMin = ::atof(params[1].c_str());
	m_eMax = ::atof(params[2].c_str());
	m_orbitalPeriod = ::atof(params[3].c_str()) * daySecs;
	m_orbitalModulation = ::atof(params[4].c_str());
	m_phi0 = ::atof(params[5].c_str());

	m_orbitalRegion.setSpectralIndex(::atof(params[6].c_str()),::atof(params[7].c_str()));
	m_orbitalRegion.setOrbitalPhase(::atof(params[8].c_str()),::atof(params[9].c_str()));

	m_diskProperties.setCycleDuration(::atof(params[10].c_str())* daySecs);
	m_diskProperties.setCycleDurationFluct(::atof(params[11].c_str()));

	m_jetProperties.setJetOnCycle(::atof(params[12].c_str()));
	m_jetProperties.setJetOnCycleFluct(::atof(params[13].c_str()));
	m_jetProperties.setJetOnDuration(::atof(params[14].c_str()));
	m_jetProperties.setJetOnDurationFluct(::atof(params[15].c_str()));



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
	double now = fmod((float)(m_currentTime+m_phi0*m_orbitalPeriod),(float)m_orbitalPeriod)/m_orbitalPeriod*2.*M_PI;
	float z = -log(1.-m_randPhase);
	m_nTurns = floor(z/(2.*M_PI*scale));
	float zp = z - 2.*M_PI*m_nTurns*scale;
	funcValue =  scale*(x - m_orbitalModulation*(sin(x+now)-sin(now))) - zp;
	derivValue = scale*(1.- m_orbitalModulation*cos(x+now));
	return;
}

double microQuasar::interval(double current_time) {
	m_currentTime = current_time - Spectrum::startTime();

	float fTime = m_currentTime;
	float jetStart = (floor(fTime/m_diskProperties.getCycleDuration()) + m_jetProperties.getJetOnCycle()) * m_diskProperties.getCycleDuration();
	float jetLength = m_jetProperties.getJetOnDuration()* m_diskProperties.getCycleDuration();
	float jetEnd = jetStart + jetLength;

	double deltaT;

	// generate times until one falls in a jet-on period. If an attempt fails, fast forward 
	// the clock to the next jet-on period and fire again.

	int i=0;
	for (i; i<100; i++) {
		m_randPhase = CLHEP::RandFlat::shoot();
		deltaT = m_orbitalPeriod/(2.*M_PI)*(rtsafe(0.,2.*M_PI,1.e-2)+2.*M_PI*m_nTurns);
		if ((m_currentTime+deltaT > jetStart) && (m_currentTime+deltaT < jetEnd)) break;
		m_currentTime = jetStart + m_diskProperties.getCycleDuration();
		jetStart += m_diskProperties.getCycleDuration();
		jetEnd += m_diskProperties.getCycleDuration();
	}
	if (i==100) std::cerr << " microQuasar::interval - exiting with max iterations " << std::endl;

	return m_currentTime - fTime + deltaT;
}


int microQuasar::OrbitalRegion::findRegion(double time, float period) {
	float timef = time;
	float phase = fmod(timef,period)/period;
	return (phase > m_minOrbitalPhase && phase < m_maxOrbitalPhase) ? 0 : 1;
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

