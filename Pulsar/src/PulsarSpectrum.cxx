#include <iostream>
#include "Pulsar/PulsarSpectrum.h"
#include "SpectObj/SpectObj.h"
#include "flux/SpectrumFactory.h"
#include "astro/JulianDate.h"

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
  m_t0 = 0.0;
  m_f0 = 1.0/m_period;
  m_f1 = 1.0/m_pdot;

  std::cout << " \n*****\nPulsarSpectrum initialized ! " << std::endl;
  std::cout << " **   Flux above 100 MeV " << m_flux << " ph/cm2/s " 
	    << " | " << m_numpeaks << " peaks " << std::endl;
  //std::cout << " **   Epoch  " << m_t0 << " , Period " << m_period << " f0 " << m_f0 << std::endl;
  //std::cout << "                         Pdot " <<  m_pdot  << " f1 " << m_f1 << std::endl; 
  std::cout << " **   Enphmin " << m_enphmin << " | Enphmax " << m_enphmax << std::endl;



  m_Pulsar   = new PulsarSim(m_flux, m_enphmin, m_enphmax, m_period, m_numpeaks);
  m_spectrum = new SpectObj(m_Pulsar->PSRPhenom(ppar1,ppar2,ppar3,ppar4),1);

  JDStartMission = new astro::JulianDate(2001,1,1,0);
  JDStartSimulation = new astro::JulianDate(2008, 2, 24, 6.5);
  m_JDCurrent = new astro::JulianDate(2008, 2, 24, 6.5);



  std::cout << " Mission Started at " << JDStartMission->seconds() <<  " ( " 
	    << JDStartMission->getGregorianDate() << " ) " << std::endl;
  std::cout << " Simulation Started at " << JDStartSimulation->seconds() <<  " ( " 
	    << JDStartSimulation->getGregorianDate() << " ) " << std::endl;
 

  //////////////////////////////////////////////////
  
}

PulsarSpectrum::~PulsarSpectrum() 
{  
  delete m_Pulsar;
  delete m_spectrum;
  delete JDStartMission;
  delete JDStartSimulation;
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
  double inte;  
  // Variable with suffix tilde are reffered to dilated system
  double timeTilde = time;
  time = timeTilde - timeTilde*m_pdot;
  //std::cout << " t0~ = " << timeTilde  << " --> t0 = " << time << std::endl;
  inte = m_spectrum->interval(time,m_enphmin);
  double nextTime = time + inte;
  double nextTimeTilde = nextTime + m_pdot*nextTime;
  //std::cout << " t1 = " << nextTime  << " --> t1~ = " << nextTimeTilde << std::endl;
  inte = nextTimeTilde - timeTilde;


  /*
  std::cout << "Next Photon at Mission Elapsed Time (sec.) " 
	    << (inte + timeTilde + JDStartSimulation->seconds() 
		- JDStartMission->seconds()) 
	    << std::endl;
  */

  return inte;
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
