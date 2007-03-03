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

  std::vector<std::string> tokens = tokenize(paramString,',');

  float daySecs = 86400.;
  
  float specOrbital1,specOrbital2=-1;
  float phaseOrbital1,phaseOrbital2=-1;

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
      m_orbitalPeriod = std::atof(token[1].c_str());
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
	curToken++;
  } 
  if(specOrbital1!=-1 && specOrbital2!=-1)
    m_orbitalRegion.setSpectralIndex(specOrbital1,specOrbital2);
  if(phaseOrbital1!=-1 && phaseOrbital2!=-1)
    m_orbitalRegion.setOrbitalPhase(phaseOrbital1,phaseOrbital2);
	/*
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
	*/

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
	
	if (m_jetStart == 0.) {
		std::pair<float,float> jet = calculateJetStart(false,fTime);
		m_jetStart = jet.first;
		m_jetEnd = jet.second;
	}


	if (fTime> m_jetEnd) {
		jet = calculateJetStart(true,fTime);
		m_jetStart = jet.first;
		m_jetEnd = jet.second;
	}

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

		double nextTime = m_currentTime+deltaT;
		if ((nextTime >= m_jetStart) && (nextTime <= m_jetEnd)) break;

		// if we're still in the current jet, don't recalculate it (because of
		// fluctuations as much as efficiency)
		if (nextTime > m_jetEnd) jet = calculateJetStart(true, nextTime);
		m_jetStart = jet.first;
		m_jetEnd = jet.second;

		m_currentTime = m_jetStart;
	}
	if (i==100) std::cerr << " microQuasar::interval - exiting with max iterations " << std::endl;

	double dT = m_currentTime - fTime + deltaT;
	return dT;
}


std::pair<double,double> microQuasar::calculateJetStart(bool nextOn, double time) {

	double jetOn = m_jetProperties.getJetOnCycle() * 
		(1. + m_jetProperties.getJetOnCycleFluct()*(0.5*CLHEP::RandFlat::shoot()-1.));

	double diskCycle = m_diskProperties.getCycleDuration() *
		(1. + m_diskProperties.getCycleDurationFluct()*(0.5*CLHEP::RandFlat::shoot()-1.));

	double jetCycle = jetOn * diskCycle;
		
	double jetLength = m_jetProperties.getJetOnDuration()* 
		(1. + m_jetProperties.getJetOnDurationFluct()*(0.5*CLHEP::RandFlat::shoot()-1.)) *
		diskCycle;

	double nJet = floor((time+jetLength+jetCycle)/diskCycle);
	if (nextOn) nJet++;
	double jetStart = jetCycle + nJet*diskCycle;
	double jetEnd = jetStart + jetLength;
	return std::make_pair(jetStart, jetEnd);
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

