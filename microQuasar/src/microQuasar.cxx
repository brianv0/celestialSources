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
#include <cctype>

#include "facilities/Util.h"

#include "flux/SpectrumFactory.h"
#include "flux/EventSource.h"

#include "microQuasar/microQuasar.h"

#include "astro/GPS.h"
#include "astro/PointingHistory.h"
#include "CLHEP/Random/RandFlat.h"

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
	float daySecs = 86400.;

	m_burstSeed = 124556789.;  // default see for burst generation

	float specOrbital1,specOrbital2=-1;
	float phaseOrbital1,phaseOrbital2=-1;

	// allow bursts to be off by default

	m_jetProperties.setJetOnDuration(1.);
	m_jetProperties.setJetOnCycle(1.);
	m_jetProperties.setJetOnCycleFluct(1.);
	m_jetProperties.setJetOnDurationFluct(1.);

	std::vector<std::string> tokens = tokenize(paramString,',');

	std::vector<std::string>::iterator curToken = tokens.begin();
	while(curToken!=tokens.end()){
		std::vector<std::string> token = tokenize(*curToken,'=');
		std::transform(token[0].begin(),token[0].end(),token[0].begin(),(int(*)(int))toupper);
		if(token[0]=="FLUX")
			m_ftot = std::atof(token[1].c_str());
		if(token[0]=="EMIN")
			m_eMin = std::atof(token[1].c_str());
		if(token[0]=="EMAX")
			m_eMax = std::atof(token[1].c_str());
		if(token[0]=="ORBITALPERIOD")
			m_orbitalPeriod = std::atof(token[1].c_str())*daySecs;
		if(token[0]=="ORBITALMODULATION")
			m_orbitalModulation = std::atof(token[1].c_str());
		if(token[0]=="ORBITALPHASE")
			m_phi0 = std::atof(token[1].c_str());
		if(token[0]=="SPECTRALORBITALREGION1")
			specOrbital1 = std::atof(token[1].c_str());
		if(token[0]=="SPECTRALORBITALREGION2")
			specOrbital2 = std::atof(token[1].c_str());
		if(token[0]=="ORBITALPHASEREGION1")
			phaseOrbital1 = std::atof(token[1].c_str());
		if(token[0]=="ORBITALPHASEREGION2")
			phaseOrbital2 = std::atof(token[1].c_str());
		if(token[0]=="DISKCYCLEDURATION")
			m_diskProperties.setCycleDuration(std::atof(token[1].c_str())*daySecs);
		if(token[0]=="DISKCYCLEFLUCTUATION")
			m_diskProperties.setCycleDurationFluct(std::atof(token[1].c_str()));
		if(token[0]=="JETONCYCLE")
			m_jetProperties.setJetOnCycle(std::atof(token[1].c_str()));
		if(token[0]=="JETONCYCLEFLUCTUATION")
			m_jetProperties.setJetOnCycleFluct(std::atof(token[1].c_str()));
		if(token[0]=="JETONDURATION")
			m_jetProperties.setJetOnDuration(std::atof(token[1].c_str()));
		if(token[0]=="JETONDURATIONFLUCTUATION")
			m_jetProperties.setJetOnDurationFluct(std::atof(token[1].c_str()));
		if(token[0]=="BURSTRANDOMSEED")
			m_burstSeed = std::atof(token[1].c_str());
		curToken++;
	} 

	// default to same spectral index in both regions
	if(specOrbital2=-1) {
		m_orbitalRegion.setSpectralIndex(specOrbital1,specOrbital1);
		m_orbitalRegion.setOrbitalPhase(0.5,1.0);
	}
	else {
		m_orbitalRegion.setSpectralIndex(specOrbital1,specOrbital2);
		m_orbitalRegion.setOrbitalPhase(phaseOrbital1, phaseOrbital2);
	}

	m_randGenBurst.setTheSeed(m_burstSeed);

	double time = 0.;
	burstPairs burstTimes;

	// get the time window for creating the burst list from the pointing history start/stop times
	astro::GPS* gps = astro::GPS::instance();
	const astro::PointingHistory& pHistory = gps->history();
    double tMin( pHistory.startTime() );
	double tMax( pHistory.endTime() );
	time = tMin;

//	std::cout << "Pointing start time (d) " << tMin/daySecs << " end time (d) " << tMax/daySecs << std::endl;

	// precompute the bursts if needed

	m_nJet = 0;
	if (m_jetProperties.getJetOnDuration() != 1.) {
		while (time < tMax) {
			burstTimes = calculateJetStart(time);
			if (burstTimes.getBurstPairs().first > tMax) break;
			else if (burstTimes.getBurstPairs().second > tMax) {
				m_jetStart = burstTimes.getBurstPairs().first;
				burstTimes = burstPairs(m_jetStart,tMax);
			}
			m_bursts.push_back(burstTimes);
			time = std::min(tMax,burstTimes.getBurstPairs().second);
//			std::cout << " Burst start (d) " << burstTimes.getBurstPairs().first/daySecs << 
//				" Burst end (d) " << time/daySecs << std::endl;
		}
	}

	m_jetStart = 0.;

	std::cerr << "microQuasar created. Total flux = " 
		<< m_ftot << " m^-2 s^-1 " << " between " 
		<< m_eMin << " MeV and "
		<< m_eMax << " MeV." << std::endl;
}

std::vector<std::string> microQuasar::tokenize(std::string params, char token){
	std::vector<std::string> tokens;
	while(params.length()!=0){
		std::string::size_type pos = params.find(token,0);
		std::string token = params.substr(0,pos);
		while(token.at(0)==' ')
			token.erase(0,1);
		while(token.at(token.length()-1)==' ')
			token.erase(token.length()-1,1);
		tokens.push_back(token);
		params.erase(0,pos);
		if(params.length()!=0)
			params.erase(0,1);
		while(params.length()!=0 && params.at(0)==' ')
			params.erase(0,1);
	}
	return tokens;
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
	// 
	double scale = m_orbitalPeriod/2./M_PI*m_ftot*EventSource::totalArea();

	// normalize 'now' so that time runs from 0-2pi, ie only look at one period
	double nowPeriodNorm = m_currentTime/m_orbitalPeriod;
	double now = (nowPeriodNorm+m_phi0)*2.*M_PI;

	double z = -log(1.-m_randPhase);
	m_nTurns = floor(z/(2.*M_PI*scale));
	double zp = z - 2.*M_PI*m_nTurns*scale;
	funcValue =  scale*(x - m_orbitalModulation*(sin(x+now)-sin(now))) - zp;
	derivValue = scale*(1.- m_orbitalModulation*cos(x+now));
	return;
}

double microQuasar::interval(double current_time) {

	m_currentTime = current_time;

	double fTime = m_currentTime;
	std::pair<double, double> jet;
	float daySecs = 86400.;


	double deltaT;
	double twoPi = 2.*M_PI;

	// generate times until one falls in a jet-on period. If an attempt fails, fast forward 
	// the clock to the next jet-on period after the projected time and fire again. Find the 
	// next jet-end period after the trial time

	int i=0;
	for (i; i<100; i++) {
		m_randPhase = CLHEP::RandFlat::shoot();
		double funcZero = rtsafe(0.,twoPi,1.e-3);
		if ( (funcZero == 0.) || (funcZero == twoPi)) {
			std::cout << "INFO: microQuasar::interval bad initial conditions to rtsafe at time " 
				<< current_time << std::endl;
			continue;
		}
		deltaT = m_orbitalPeriod/twoPi*(funcZero + twoPi*m_nTurns);

		// if steady source, don't worry about artificial jet-on boundaries
		if (m_jetProperties.getJetOnDuration() == 1.) break;
/*		Handling 4 conditions for outbursts:
			1. current time is before first outburst: move clock to first outburst and get new time
			2. during an outburst - accept time
			3. between outburst - same as 1
			4. after last outburst - terminate
*/
//		jet = getJetStart(fTime);
//		m_jetStart = jet.first;
//		m_jetEnd = jet.second;

		double nextTime = m_currentTime+deltaT;
		if (inJet(nextTime)) break;

		std::vector<burstPairs>::iterator jet = getJetStart(nextTime);
		if (jet != m_bursts.end())	m_currentTime = (*jet).getBurstPairs().first;
		else {
			deltaT = 1.E9;
			break;
		}
	}
	if (i==100) std::cerr << " microQuasar::interval - exiting with max iterations " << std::endl;

	double dT = m_currentTime - fTime + deltaT;
//	std::cout << "input t (d) " << fTime/daySecs << " interval (sec) " << dT << std::endl;
	return dT;
}

bool microQuasar::inJet(double time) 
{
    std::vector<burstPairs>::iterator burstIter = std::find_if(m_bursts.begin(), m_bursts.end(), InBurst(time));
    if (burstIter!= m_bursts.end()) return true;
	return false;
}

std::vector<microQuasar::burstPairs>::iterator microQuasar::getJetStart(double time) 
{
    std::vector<burstPairs>::iterator burstIter = std::find_if(m_bursts.begin(), m_bursts.end(), NextBurst(time));
    return burstIter;
}
microQuasar::burstPairs microQuasar::calculateJetStart(double time) {

	double randJet = 0.5*(2.*m_randGenBurst.flat()-1.);
	double jetOn = m_jetProperties.getJetOnCycle() * 
		(1. + m_jetProperties.getJetOnCycleFluct()* randJet);

	randJet = 0.5*(2.*m_randGenBurst.flat()-1.);
	double diskCycle = m_diskProperties.getCycleDuration() * 
		(1. + m_diskProperties.getCycleDurationFluct()* randJet);

	double jetCycle = jetOn * diskCycle;

	randJet = 0.5*(2.*m_randGenBurst.flat()-1.);
	double jetLength = m_jetProperties.getJetOnDuration()* 
		(1. + m_jetProperties.getJetOnDurationFluct()* randJet) *
		diskCycle;

//	if (m_nJet== 0) m_nJet = floor((time+jetLength+jetCycle)/diskCycle);
	if (m_nJet== 0) {
		m_nJet = time/diskCycle;
		m_cycleStart = time;
	}
	else {
		m_nJet++;
		m_cycleStart += diskCycle;
	}

	double jetStart = jetCycle + m_cycleStart;  
	double jetEnd = jetStart + jetLength;
	return burstPairs(jetStart,jetEnd);
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
