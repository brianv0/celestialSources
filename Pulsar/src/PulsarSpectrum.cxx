#include <iostream>
#include "Pulsar/PulsarSpectrum.h"
#include "Pulsar/PulsarConstants.h"
#include "SpectObj/SpectObj.h"
#include "flux/SpectrumFactory.h"
#include "astro/JulianDate.h"
#include <cmath>
#include <fstream>

using namespace cst;

ISpectrumFactory &PulsarSpectrumFactory() 
 {
   static SpectrumFactory<PulsarSpectrum> myFactory;
   return myFactory;
 }

PulsarSpectrum::PulsarSpectrum(const std::string& params)
  : m_params(params)
{

  // Set to 0 all variables;
  m_RA = 0.0;
  m_dec = 0.0;
  m_l = 0.0;
  m_b = 0.0;
  m_flux = 0.0;
  m_period = 0.0;
  m_pdot = 0.0;
  m_t0 = 0.0;
  m_phi0 = 0.0;
  m_model = 0;
  m_seed = 0;

  char* pulsar_root = ::getenv("PULSARROOT");
  char sourceFileName[80];
  std::ifstream PulsarDataTXT;

  //Read from XML file

  m_PSRname    = parseParamList(params,0).c_str();
  m_model      = std::atoi(parseParamList(params,1).c_str());

  m_enphmin    = std::atof(parseParamList(params,2).c_str());
  m_enphmax    = std::atof(parseParamList(params,3).c_str());

  m_seed       = std::atoi(parseParamList(params,4).c_str()); //Seed for random
  m_numpeaks   = std::atoi(parseParamList(params,5).c_str()); //Number of peaks

  double ppar1 = std::atof(parseParamList(params,6).c_str());
  double ppar2 = std::atof(parseParamList(params,7).c_str());
  double ppar3 = std::atof(parseParamList(params,8).c_str());
  double ppar4 = std::atof(parseParamList(params,9).c_str());

  //Read from PulsarDataList.txt
  
  sprintf(sourceFileName,"%s/data/PulsarDataList.txt",pulsar_root);  
  std::cout << "\nOpening Pulsar Datalist File : " << sourceFileName << std::endl;
  PulsarDataTXT.open(sourceFileName, std::ios::in);
  
  if (! PulsarDataTXT.is_open()) 
    {
      std::cout << "Error opening Datalist file " << sourceFileName  
		<< " (check whether $PULSARROOT is set" << std::endl; 
      exit (1);
    }
  
  char aLine[200];  
  PulsarDataTXT.getline(aLine,200); 

  char tempName[15] = "";
  
  while ((std::string(tempName) != m_PSRname) && ( PulsarDataTXT.eof() != 1))
    {
      PulsarDataTXT >> tempName >> m_RA >> m_dec >> m_l >> m_b >> m_flux >> m_period >> m_pdot >> m_t0 >> m_phi0;
      //  std::cout << tempName  << m_RA << m_dec << m_l << m_b << m_flux << m_period << m_pdot << m_t0 << m_phi0 << std::endl;  
    } 
  
  if (std::string(tempName) == m_PSRname)
    {
      std::cout << "Pulsar " << m_PSRname << " found in Datalist file! " << std::endl;
    }
  else
    {
      std::cout << "ERROR! Pulsar " << m_PSRname << " NOT found in Datalist.File...Exit. " << std::endl;
      exit(1);
    }
  
  m_f0 = 1.0/m_period;
  m_f1 = -m_pdot/(m_period*m_period);
  
 
  //Writing out Pulsar Info
  std::cout << " \n********   PulsarSpectrum initialized !   ********" << std::endl;
  std::cout << "**   Name : " << m_PSRname << std::endl;
  std::cout << "**   Position : (RA,Dec)=(" << m_RA << "," << m_dec 
	    << ") ; (l,b)=(" << m_l << "," << m_b << ")" << std::endl; 
  std::cout << "**   Flux above 100 MeV : " << m_flux << " ph/cm2/s " << std::endl;
  std::cout << "**   Number of peaks : " << m_numpeaks << std::endl;
  std::cout << "**   Epoch (MJD) :  " << m_t0 << std::endl;
  std::cout << "**   Phi0 (at Epoch t0) : " << m_phi0 << std::endl;
  std::cout << "**   Period : " << m_period << " s. | f0: " << m_f0 << std::endl;
  std::cout << "**   Pdot : " <<  m_pdot  << " | f1: " << m_f1 << std::endl; 
  std::cout << "**   Enphmin : " << m_enphmin << " keV | Enphmax: " << m_enphmax << " keV" << std::endl;
  std::cout << "**   Mission started at (MJD) : " << StartMissionDateMJD << " (" 
	    << (StartMissionDateMJD+JDminusMJD)*SecsOneDay << " sec.) - Jan,1 2001 00:00.00" << std::endl;
  std::cout << "**************************************************" << std::endl;


  m_JDCurrent = new astro::JulianDate(2001, 1, 1, 0.0);
  m_Pulsar    = new PulsarSim(m_PSRname, m_seed, m_flux, m_enphmin, m_enphmax, m_period, m_numpeaks);

  if (m_model == 1)
    {
      std::cout << "\n**   Model chosen : " << m_model << " --> Using Phenomenological Pulsar Model " << std::endl;  
      m_spectrum = new SpectObj(m_Pulsar->PSRPhenom(ppar1,ppar2,ppar3,ppar4),1);
      m_spectrum->SetAreaDetector(EventSource::totalArea());
      std::cout << "**  Effective Area set to : " << m_spectrum->GetAreaDetector() << " m2 " << std::endl; 
    }
  else 
    {
      std::cout << "ERROR!  Model choice not implemented " << std::endl;
      exit(1);
    }

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
  if ((int(time) % 10000) < 10)
    std::cout << "Time reached is: " << time << " seconds from Start " << std::endl;


  if (m_pdot == 0.0) // Case of periodic pulsar
    {
      time = timeTilde + m_phi0*m_period;
    }
  else
    {
      time = (m_period/m_pdot)*log(1+(m_pdot/m_period)*timeTilde) + m_phi0*m_period; 
    }
  
  double nextTime = time + m_spectrum->interval((time - m_period*int(time/m_period)),m_enphmin);

  double nextTimeTilde = 0.0;

  if (m_pdot == 0.0) //Case of periodic pulsar
    {
      nextTimeTilde = nextTime - m_phi0*m_period;
    }
  else
    {
      nextTimeTilde = (m_period/m_pdot)*(exp((m_pdot/m_period)*(nextTime-m_phi0*m_period))-1);
    }

  double ph =  m_phi0 + m_f0*(timeTilde) +0.5*m_f1*(timeTilde)*(timeTilde);
  double intph = 0;

  
  if ((fabs(modf(ph,&intph) - ((time -m_period*int(time/m_period))/m_period)) > 2.5e-3)
      && (timeTilde - (StartMissionDateMJD - m_t0)*SecsOneDay > 0.0))
    {
      
      std::cout << "\n****t0~ = " << timeTilde  << " (RefStart="
		<< timeTilde - (StartMissionDateMJD - m_t0)*SecsOneDay  << ")" << std::endl;
      std::cout << " -->t0  = " << time << " (Fract = " << (time  -m_period*int(time/m_period)) << " ), phase " 
		<<  (time - m_period*int(time/m_period))/m_period << std::endl;
      std::cout << "\n****t1~ = " << nextTimeTilde  << " (RefStart="
		<< nextTimeTilde - (StartMissionDateMJD - m_t0)*SecsOneDay  << ")" << std::endl;
      std::cout << " -->t1  = " << nextTime << " (Fract = " << (nextTime -m_period*int(nextTime/m_period)) << " ), phase " 
		<< (nextTime -m_period*int(nextTime/m_period))/m_period << std::endl;

      
      std::cout << "\n\n** " <<  m_PSRname << " Warning: 2 Order Approx (LAT-ST pulsePhase):  ph0 " 
		<< ph <<  " " << modf(ph,&intph) << std::endl;
      std::cout << "**          Phase difference : ph0 " 
		<< fabs(modf(ph,&intph) - ((time -m_period*int(time/m_period))/m_period)) << std::endl;

      double ph =  m_phi0 + m_f0*(nextTimeTilde) +0.5*m_f1*(nextTimeTilde)*(nextTimeTilde);
      std::cout << "**          Phase difference : ph1 " 
		<< fabs(modf(ph,&intph) - ((nextTime -m_period*int(nextTime/m_period))/m_period)) << std::endl;


    }


  return nextTimeTilde - timeTilde;
}

double PulsarSpectrum::energy(double time)
{
  return m_spectrum->energy(time,m_enphmin)*1.0e-3; //MeV
}

std::string PulsarSpectrum::parseParamList(std::string input, int index)
{
  std::vector<std::string> output;
  unsigned int i=0;
  
  for(;!input.empty() && i!=std::string::npos;){
   
    i=input.find_first_of(",");
    std::string f = ( input.substr(0,i).c_str() );
    output.push_back(f);
    input= input.substr(i+1);
  } 

  if(index>=output.size()) return "";
  return output[index];
};



