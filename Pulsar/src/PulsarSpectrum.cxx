#include <iostream>
#include "Pulsar/PulsarSpectrum.h"
#include "Pulsar/PulsarConstants.h"
#include "SpectObj/SpectObj.h"
#include "flux/SpectrumFactory.h"
#include "astro/JulianDate.h"
#include <cmath>

using namespace cst;

ISpectrumFactory &PulsarSpectrumFactory() 
 {
   static SpectrumFactory<PulsarSpectrum> myFactory;
   return myFactory;
 }

PulsarSpectrum::PulsarSpectrum(const std::string& params)
  : m_params(params)
{
  std::cout<<m_flux<<std::endl;
  m_flux = parseParamList(params,1); //Flux above 100 MeV in ph/cm2/s
  m_period  = parseParamList(params,2); //Period in sec.
  m_pdot = parseParamList(params,3); //Pdot
  m_numpeaks = parseParamList(params,4); //Number of peaks

  double ppar1 = parseParamList(params,5);
  double ppar2 = parseParamList(params,6);
  double ppar3 = parseParamList(params,7);
  double ppar4 = parseParamList(params,8);
  m_enphmin = parseParamList(params,9);
  m_enphmax = parseParamList(params,10);
  m_phi0 = 0.032;
  m_t0 = 51900.0;
  m_f0 = 1.0/m_period;
  m_f1 = -m_pdot/(m_period*m_period);
  
  std::cout << " \n********   PulsarSpectrum initialized !   ********" << std::endl;
  std::cout << "**   Flux above 100 MeV : " << m_flux << " ph/cm2/s " << std::endl;
  std::cout << "**   # of peaks : " << m_numpeaks << std::endl;
  std::cout << "**   Epoch (MJD) :  " << m_t0 << std::endl;
  std::cout << "**   Phi0 (at Epoch t0) : " << m_phi0 << std::endl;
  std::cout << "**   Period : " << m_period << " s. | f0: " << m_f0 << std::endl;
  std::cout << "**   Pdot : " <<  m_pdot  << " | f1: " << m_f1 << std::endl; 
  std::cout << "**   Enphmin : " << m_enphmin << " keV | Enphmax: " << m_enphmax << " keV" << std::endl;
  std::cout << "**   Mission started at (MJD) : " << StartMissionDateMJD << " (" 
	    << (StartMissionDateMJD+JDminusMJD)*SecsOneDay << " sec.) - Jan,1 2001 00:00.00" << std::endl;
  std::cout << "**************************************************" << std::endl;

  m_Pulsar   = new PulsarSim(m_flux, m_enphmin, m_enphmax, m_period, m_numpeaks);
  m_spectrum = new SpectObj(m_Pulsar->PSRPhenom(ppar1,ppar2,ppar3,ppar4),1);
  m_JDCurrent = new astro::JulianDate(2008, 2, 24, 6.5);

  //////////////////////////////////////////////////
  
}

PulsarSpectrum::~PulsarSpectrum() 
{  
  delete m_Pulsar;
  delete m_spectrum;
  delete m_JDCurrent;
}

//return flux, given a time
double PulsarSpectrum::flux(double time) const
{
  double flux;	  
  flux = m_spectrum->flux(time,m_enphmin);
  return flux;
}

double PulsarSpectrum::interval(double time)
{  
  double timeTilde = time + (StartMissionDateMJD - m_t0)*SecsOneDay;
  time = (m_period/m_pdot)*log(1+(m_pdot/m_period)*timeTilde) + m_phi0*m_period; 
  double nextTime = time + m_spectrum->interval((time - m_period*int(time/m_period)),m_enphmin);
  double nextTimeTilde = (m_period/m_pdot)*(exp((m_pdot/m_period)*(nextTime-m_phi0*m_period))-1);



  std::cout << "\n****t0~ = " << timeTilde  << " (RefStart="
	    << timeTilde - (StartMissionDateMJD - m_t0)*SecsOneDay  << ")" << std::endl;
  std::cout << " -->t0  = " << time << " (Fract = " << (time - m_period*int(time/m_period)) << " ), phase " 
	    <<  (time - m_period*int(time/m_period))/m_period << std::endl;
  std::cout << " ----> Interval " << inte << std::endl;
  std::cout << "\n****t1~ = " << nextTimeTilde  << " (RefStart="
	    << nextTimeTilde - (StartMissionDateMJD - m_t0)*SecsOneDay  << ")" << std::endl;
  std::cout << " -->t1  = " << nextTime << " (Fract = " << (nextTime - m_period*int(nextTime/m_period)) << " ), phase " 
	    << (nextTime - m_period*int(nextTime/m_period))/m_period << std::endl;
 
  double ph = ph = m_phi0 + m_f0*(timeTilde) +0.5*m_f1*(timeTilde)*(timeTilde);
  double intph = 0;



  std::cout << "\n\n *** theory : ph0 " << ph <<  " " << modf(ph,&intph) << std::endl;
  return nextTimeTilde - timeTilde;
}

double PulsarSpectrum::energy(double time)
{
  return m_spectrum->energy(time,m_enphmin)*1.0e-3; //MeV
}

double PulsarSpectrum::parseParamList(std::string input, int index)
{
  std::vector<double> output;
  unsigned int i=0;
  for(;!input.empty() && i!=std::string::npos;){
    double f = ::atof( input.c_str() );
    output.push_back(f);
    i=input.find_first_of(",");
    input= input.substr(i+1);
  } 
  if(index>=output.size()) return 0.0;
  return output[index];
};
