#include <iostream>
#include "Pulsar/PulsarSpectrum.h"
#include "Pulsar/PulsarConstants.h"
#include "SpectObj/SpectObj.h"
#include "flux/SpectrumFactory.h"
#include "astro/JulianDate.h"
#include <cmath>
#include <fstream>
#include <iomanip>


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
  m_p2dot = 3.0e-18;
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
  m_f2 = 2*pow((m_pdot/m_period),2.0)/m_period - m_p2dot/(m_period*m_period);
  
 
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
  double timeTilde = time + (StartMissionDateMJD)*SecsOneDay;
  if ((int(timeTilde - (StartMissionDateMJD)*SecsOneDay) % 20000) < 2.0)
    std::cout << "**  Time reached is: " << timeTilde-(StartMissionDateMJD)*SecsOneDay
	      << " seconds from Mission Start " << std::endl;
  
  double initTurns = getTurns(timeTilde); //Turns made at this time 
  double intPart; //Integer part
  double tStart = modf(initTurns,&intPart)*m_period; // Start time for interval
  double inte = m_spectrum->interval(tStart,m_enphmin); //deltaT (in system where Pdot = 0
  double inteTurns = inte/m_period; // conversion to # of turns
  double totTurns = initTurns + inteTurns; //Total turns at the nextTimeTilde
  //std::cout <<  std::setprecision(30) << "**  Inte " << inte 
  //    << " initTurns " << initTurns << " totTurns " << totTurns << std::endl;


  double finPh =  modf(initTurns,&intPart);
  finPh = modf(finPh + inteTurns,&intPart); //finPh is the expected phase corresponding to the nextTimeTilde
  double nextTimeTilde = retrieveNextTimeTilde(timeTilde, totTurns, 1e-4);
  //std::cout << std::setprecision(30) << "phCalc " << modf(getTurns(nextTimeTilde),&intPart) << " exp " << finPh << std::endl;
  //std::cout << std::setprecision(30) << "       diff " << fabs(modf(getTurns(nextTimeTilde),&intPart)-finPh) << std::endl;


  /* ///Cut the old code...

   if (m_pdot != 0.0) 
    {
      // The exact formula is:
      //     time = (m_period/m_pdot)*log(1.0+(m_pdot/m_period)*(timeTilde)) + m_phi0*m_period;
      // but we use an approximation in order to avoid errors due to finite precision in log (1+e)
      // where e is of order of 10^-12. double has around 20 decimal positions.We use instead the 
      // expansion of log(1+x)
       time = m_phi0*m_period;
       for (int n=1; n<4; n++)
	{
	  time = time + ((m_period/m_pdot)*pow(-1,n+1)*pow(((m_pdot/m_period)*timeTilde),n))/n;
	}
      std::cout << std::setprecision(20) << "time with series " << time << std::endl;
      double time1 = timeTilde -(m_pdot/(2*m_period))*timeTilde*timeTilde + (m_pdot*m_pdot)/(3*m_period*m_period)*(timeTilde*timeTilde*timeTilde);
       time = time1;
       std::cout << std::setprecision(20) << "time with 2 order " << time1 << std::endl;
      std::cout << std::setprecision(20) << "diff " << time1-time << std::endl;


    } 
  else // Case of exact periodic pulsar (pdot = 0)
    {
      time = timeTilde + m_phi0*m_period;
    }
  //Get interval in the reference system where pdot = 0;
   double nextTime = time + m_spectrum->interval((time - m_period*int(time/m_period)),m_enphmin);
  double nextTimeTilde = 0.0;

   if (m_pdot != 0.0) //Case of periodic pulsar
    {
      // The exact formula is:
      //     nextTimeTilde = (m_period/m_pdot)*(exp((m_pdot/m_period)*(nextTime-m_phi0*m_period))-1.0);
      // but we use an approximation in order to avoid errors due to finite precision in log (1+e)
      // where e is of order of 10^-12. double has around 20 decimal positions.We use instead the 
      // expansion of exp(x)
 
      
      for (int n=1; n<3;n++)
	{
	  double fact = 1.0;
	  for (int i = 0; i <n; i++) //compute factorial
	    {
	      if (i == 0) 
 		fact = 1;
 	      else
 		fact = fact*(i);
	    }
 	  nextTimeTilde = nextTimeTilde + (m_period/m_pdot)*(pow(((m_pdot/m_period)*(nextTime-m_phi0*m_period)),n)/fact);
	}
    }
  else // case of exactly periodic pulsar (pdot =0)
    {
      nextTimeTilde = nextTime - m_phi0*m_period;
    }
  
  //This is the formula for the phase computation in the Pulsar Analysys tool. Here for comparison.
  double ph0 =  m_phi0 + m_f0*(timeTilde) +0.5*m_f1*(timeTilde)*(timeTilde);
   double ph1 =  m_phi0 + m_f0*(nextTimeTilde) +0.5*m_f1*(nextTimeTilde)*(nextTimeTilde);
  double intph0 = 0;
  double intph1 = 0;

  //The precision in time is set to usec
    
  //if ((timeTilde - (StartMissionDateMJD - m_t0)*SecsOneDay > 0.0)  //check if Dt is > of 10us
  //  && ((fabs(modf(ph0,&intph0) - ((time -m_period*int(time/m_period))/m_period))*m_period) > 5.0e-5))
	  //  || ((fabs(modf(ph1,&intph1) - ((nextTime -m_period*int(nextTime/m_period))/m_period))*m_period) > 5.0e-5)))
  //  {
  
      std::cout << "\n\n** " <<  m_PSRname << " Warning: 2 Order Approx (LAT-ST pulsePhase): " <<std::endl; 
      std::cout << std::setprecision(20) << "\n**** t0~ = " << timeTilde  << " (RefStart="
		<< timeTilde - (StartMissionDateMJD - m_t0)*SecsOneDay  << ")" << std::endl;
      std::cout << std::setprecision(20) << "  -->t0 = " << time << " (Fract = " << (time  -m_period*int(time/m_period)) << " ), phase " 
		<<  (time - m_period*int(time/m_period))/m_period << std::endl;
            
      std::cout << std::setprecision(20)<< "\n  ==>t1 = " << nextTime << " (Fract = " << (nextTime -m_period*int(nextTime/m_period)) << " ), phase " 
		<< (nextTime -m_period*int(nextTime/m_period))/m_period << std::endl;
      std::cout << std::setprecision(20) << "  -->t1~ = " << nextTimeTilde  << " (RefStart="
		<< nextTimeTilde - (StartMissionDateMJD - m_t0)*SecsOneDay  << ")" << std::endl;
      
      std::cout <<  std::setprecision(16)   << "**  ph0(th) : " << modf(ph0,&intph0) << " ph0(sim) : " 
		<< ((time -m_period*int(time/m_period))/m_period)
		<< " deltaPh0 " << fabs(modf(ph0,&intph0) - ((time -m_period*int(time/m_period))/m_period)) << std::endl;
      std::cout	<<  std::setprecision(16)  <<  "     deltat0 " 
		<<  m_period*fabs(modf(ph0,&intph0) - ((time -m_period*int(time/m_period))/m_period)) << std::endl;
  
      std::cout <<  std::setprecision(16)   << "**  ph1(th) : " << modf(ph1,&intph1) << " ph1(sim) : " 
		<< ((nextTime -m_period*int(nextTime/m_period))/m_period)
		<< " deltaPh1 " << fabs(modf(ph1,&intph1) - ((nextTime -m_period*int(nextTime/m_period))/m_period))  << std::endl;
      std::cout	<<  std::setprecision(16) <<   "     deltat1 " 
		<<  m_period*fabs(modf(ph1,&intph1) - ((nextTime -m_period*int(nextTime/m_period))/m_period)) << " \n\n" << std::endl;
      
      //}


  */ // end of Cut...


   return nextTimeTilde - timeTilde;
}

/////////////////////////////////////////////
double PulsarSpectrum::getTurns( double time )
{
  double dt = time - m_t0*SecsOneDay;
  return m_phi0 + m_f0*dt + 0.5*m_f1*dt*dt + ((m_f2*dt*dt*dt)/6.0);
}

////////////////////////////////////////////
double PulsarSpectrum::retrieveNextTimeTilde( double tTilde, double totalTurns, double err )
{
  //err = 1e-6;
  double tStep = m_period;
  double tempTurns = 0.0;
  double tTildeHigh = 0.0;
  double tTildeLow = tTilde;
  //  std::cout << " Retrieving Period " << std::endl;
  //  std::cout << std::setprecision(30) << " tTilde " << tTilde << " totalTurns " << totalTurns << " err " << err << std::endl;
  //First part to find Extrema...

  while (tempTurns < totalTurns)
    {
      tTilde = tTilde + tStep;
      tempTurns = getTurns(tTilde);
    }
  tTildeHigh = tTilde;
  //std::cout << "\n **First Iter: Low is " << tTildeLow << " High is " << tTildeHigh << " -->turns " << tempTurns <<std::endl;


  //Iterative procedure to find the correct nextTime

 tStep = (tTildeHigh-tTildeLow)/2;

  while (fabs(tempTurns - totalTurns ) > err)
    {

      tTilde = tTildeLow+tStep;
      tempTurns = getTurns(tTilde);
      if ((tempTurns-totalTurns) > 0.0)
	{
	  tTildeHigh = tTilde;
	}
      else 
	{
	  tTildeLow = tTilde;
	}
      //   std::cout <<  std::setprecision(30) << " Low is " << tTildeLow << " High is " << tTildeHigh << " -->turns " << tempTurns <<std::endl;;
      tStep = tStep/2;
    }

  return tTilde;
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



