//////////////////////////////////////////////////
// File PulsarSim.cxx
// contains the code for the implementation of the models
//////////////////////////////////////////////////

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

/////////////////////////////////////////////////
ISpectrumFactory &PulsarSpectrumFactory() 
 {
   static SpectrumFactory<PulsarSpectrum> myFactory;
   return myFactory;
 }

/////////////////////////////////////////////////
/*!
 * \param params string containing the XML parameters
 *
 * This method takes from the XML file the name of the pulsar to be simulated, the model
 * to be used and the parameter specific for the model choosen. Then extract from the
 * txt DataList file the characteristics of the pulsar choosen (period, flux, period derivatives
 * ,ephemerides, flux, etc...)
 * Then it calls the PulsarSim class in order to have the TH2D histogram.
 * For more informations and fot a brief tutorial about the use of PulsarSpectrum  please see:
 * <br>
 * <a href="#dejager02">http://www.pi.infn.it/~razzano/Pulsar/PulsarSpectrumTutorial/PsrSpectrumTut.html </a>
 */
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
  m_p2dot = 0.0;
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
      PulsarDataTXT >> tempName >> m_RA >> m_dec >> m_l >> m_b >> m_flux >> m_period >> m_pdot >> m_p2dot >> m_t0 >> m_phi0;
      //  std::cout << tempName  << m_RA << m_dec << m_l << m_b << m_flux 
      //  << m_period << m_pdot << m_t0 << m_phi0 << std::endl;  
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
  
 
  //Writes out Pulsar Info
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
  std::cout << "**   P2dot : " <<  m_p2dot  << " | f2: " << m_f2 << std::endl; 
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
  
}

/////////////////////////////////////////////////
PulsarSpectrum::~PulsarSpectrum() 
{  
  delete m_Pulsar;
  delete m_spectrum;
  delete m_JDCurrent;
}

/////////////////////////////////////////////////
double PulsarSpectrum::flux(double time) const
{
  double flux;	  
  flux = m_spectrum->flux(time,m_enphmin);
  return flux;
}


/////////////////////////////////////////////////
/*!
 * \param time time to start the computing of the next photon
 *
 * This method find the time when the next photon will come. The decorrections for the ephemerides
 * are taken into account. The parameters used for this computation are:
 * <ul>
 * <li> <i>Epoch</i>, espressed in MJD;</li>
 * <li> <i>Initial Phase Phi0</i>;</li>
 * <li> <i>First Period derivative</i>;</li>
 * <li> <i>Second Period derivative</i>;</li>
 * </ul>
 * <br>
 * The interval method is taken from SpectObj. 
 */
double PulsarSpectrum::interval(double time)
{  
  double timeTilde = time + (StartMissionDateMJD)*SecsOneDay;
  if ((int(timeTilde - (StartMissionDateMJD)*SecsOneDay) % 20000) < 1.5)
    std::cout << "**  Time reached is: " << timeTilde-(StartMissionDateMJD)*SecsOneDay
	      << " seconds from Mission Start " << std::endl;
  double initTurns = getTurns(timeTilde); //Turns made at this time 
  double intPart; //Integer part
  double tStart = modf(initTurns,&intPart)*m_period; // Start time for interval
  //determines the interval to the next photon (at nextTime) in the case of period constant.
  double inte = m_spectrum->interval(tStart,m_enphmin); //deltaT (in system where Pdot = 0
  double inteTurns = inte/m_period; // conversion to # of turns
  double totTurns = initTurns + inteTurns; //Total turns at the nextTimeTilde

  double finPh =  modf(initTurns,&intPart);
  finPh = modf(finPh + inteTurns,&intPart); //finPh is the expected phase corresponding to the nextTimeTilde
  double nextTimeTilde = retrieveNextTimeTilde(timeTilde, totTurns, (1e-6/m_period));
  
  //std::cout << std::setprecision(30) << "phCalc " << modf(getTurns(nextTimeTilde),&intPart) << " exp " << finPh << std::endl;
  //std::cout << std::setprecision(30) << "       diff " << fabs(modf(getTurns(nextTimeTilde),&intPart)-finPh) << std::endl;

   return nextTimeTilde - timeTilde;
}

/////////////////////////////////////////////
/*!
 * \param time time where the number of turns is computed
 *
 * This method compute the number of turns made by the pulsar according to the ephemerides.
 * <br>
 * <blockquote>
 * f(t) = f<sub>0</sub> + f<sub>1</sub>(t - t<sub>0</sub>) +
 *      f<sub>2</sub>/2(t - t<sub>0</sub>)<sup>2</sup>.
 * <br>The number of turns is related to the ephemerides as:<br>
 * dN(t) = f<sub>0</sub> + f<sub>1</sub>(t-t<sub>0</sub>) + 1/2f<sub>2</sub>(t-t<sub>0</sub>)<sup>2</sup>
 * </blockquote>
 * <br>where 
 * <blockquote>
 *<ul>
 * <li> t<sub>0</sub> is an ephemeris epoch.In this simulator epoch has to be expressed in MJD</li>
 * <li>f<sub>0</sub> pulse frequency at time t<sub>0</sub> (usually given in Hz)</li>,
 * <li>f<sub>1</sub> the 1st time derivative of frequency at time t< sub>0</sub> (Hz s<sup>-1</sup>);</li>
 * <li>f<sub>2</sub> the 2nd time derivative of frequency at time t<sub>0</sub> (Hz s<sup>-2</sup>).</li>
 * </ul>
 * </blockquote>
 */
double PulsarSpectrum::getTurns( double time )
{
  double dt = time - m_t0*SecsOneDay;
  return m_phi0 + m_f0*dt + 0.5*m_f1*dt*dt + ((m_f2*dt*dt*dt)/6.0);
}

/////////////////////////////////////////////
/*!
 * \param tTilde Time from where retrieve the nextTime;
 * \param totalTurns Number of turns completed by the pulsar at nextTime;
 * \param err Tolerance between totalTurns and the number of turns at nextTime
 *
 * <br>In this method a recursive way is used to find the <i>nextTime</i>. Starting at <i>tTilde</i>, the method returns 
 * nextTime, the time where the number of turns completed by the pulsar is equal to totalTurns.  
 */
double PulsarSpectrum::retrieveNextTimeTilde( double tTilde, double totalTurns, double err )
{
  //err = 1e-6;
  double tStep = m_period;
  double tempTurns = 0.0;
  double tTildeHigh = 0.0;
  double tTildeLow = tTilde;
 
  //std::cout << std::setprecision(30) << " tTilde " << tTilde << " totalTurns " << totalTurns << " err " << err << std::endl;

  //First part to find Extrema...
  while (tempTurns < totalTurns)
    {
      tTilde = tTilde + tStep;
      tempTurns = getTurns(tTilde);
    }
  tTildeHigh = tTilde;
  
  //Iterative procedure to find the correct nextTime
  tStep = (tTildeHigh-tTildeLow)/2;
  
  int nIter = 0;
  
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
      // std::cout <<  std::setprecision(30) << " Low is " << tTildeLow << " High is " << tTildeHigh 
      //           << " -->turns " << tempTurns <<std::endl;;
      tStep = tStep/2;
      nIter++;

      //If the procedure not converge within the err,  the tolerance is amplyfied.
      if (nIter == 50)
	{
	  //	  	  std::cout << std::setprecision(30) << " Warning!! Amplifying tolerance for convergence at time " 
	  //		    << tTilde - (StartMissionDateMJD)*SecsOneDay << std::endl;
	  //	  std::cout <<  std::setprecision(30) << " Low is " << tTildeLow << " High is " << tTildeHigh 
	  //    << " -->turns " << tempTurns << "(total " << totalTurns << " ) " << std::endl;
	  //std::cout << "    diff " << fabs(tempTurns-totalTurns) << " , err " << err << std::endl;
	  //std::cout << " Error from " << err ; 
	  err = err*5.0;
	  //std::cout << " to " << err << " diff " << fabs(tempTurns-totalTurns) <<std::endl;
	  nIter = 0;
	}
    }

  return tTilde;
}

/////////////////////////////////////////////
double PulsarSpectrum::energy(double time)
{
  return m_spectrum->energy(time,m_enphmin)*1.0e-3; //MeV
}

/////////////////////////////////////////////
/*!
 * \param input String to be parsed;
 * \param index Position of the parameter to be extracted from input;
 *
 * <br>
 * From a string contained in the XML file a parameter is extracted according to the position 
 * specified in <i>index</i> 
 */
std::string PulsarSpectrum::parseParamList(std::string input, unsigned int index)
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
