/////////////////////////////////////////////////
// File PulsarSpectrum.cxx
// Implementation of PulsarSpectrum class
//////////////////////////////////////////////////

#include "Pulsar/PulsarSpectrum.h"
#include "Pulsar/PulsarConstants.h"
#include "SpectObj/SpectObj.h"
#include "flux/SpectrumFactory.h"

#define DEBUG 0

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
 * to be used and the parameters specific for the choosen model. Then extracts from the
 * TXT DataList file the characteristics of the choosen pulsar (e.g. period, flux, period derivatives
 * ,ephemerides, flux, etc...)
 * The names of DataList files are defined in the file <i>PulsarDataList.txt</i><br>
 * The parameters are:<br>
 * <ul>
 *  <li> Pulsar name;
 *  <li> Flux, expressed in ph/cm2/s (E>100MeV);
 *  <li> Ephemerides type (P=period and derivatives;F=frequency and derivatives)
 *  <li> Period,pdot,p2dot or f0,f1,f2;
 *  <li> Start of ephemerides validity range (in MJD);
 *  <li> Reference Epoch, t0 (in MJD);
 *  <li> End of ephemerides validity range (in MJD);
 * </ul> 
 * The parameters taken from XML file are:<br>
 * <ul>
 *  <li> Pulsar name;
 *  <li> Position (RA,dec) in degrees;
 *  <li> Minimum and maximum energy of extracted photons (in keV);
 *  <li> Model choosen (default = 1, phenomenological);
 *  <li> Random seed;
 *  <li> Model-dependend parameters;
 * </ul>
 * Then it calls the PulsarSim class in order to have the TH2D histogram.
 * For more informations and for a brief tutorial about the use of PulsarSpectrum you can look at:
 * <br>
 * <a href="#dejager02">http://www.pi.infn.it/~razzano/Pulsar/PulsarSpTutor/PulsarSpTutor.htm </a>
 */


PulsarSpectrum::PulsarSpectrum(const std::string& params)
  : m_solSys(astro::SolarSystem::EARTH), m_params(params) 
{

  // Reset all variables;
  m_RA = 0.0;
  m_dec = 0.0;
  m_l = 0.0;
  m_b = 0.0;
  m_flux = 0.0;
  m_ephemType = "";
  m_period = 0.0;
  m_pdot = 0.0;
  m_p2dot = 0.0;
  m_f0 = 0.0;
  m_f1 = 0.0;
  m_f2 = 0.0;
  m_N0 = 0.0;
  m_t0Init = 0.0;
  m_t0 = 0.0;
  m_t0End = 0.0;  
  m_phi0 = 0.0;
  m_model = 0;
  m_seed = 0;

  //Read from XML file
  m_PSRname    = parseParamList(params,0).c_str();            // Pulsar name
  m_RA = std::atof(parseParamList(params,1).c_str());         // Pulsar Right Ascension
  m_dec = std::atof(parseParamList(params,2).c_str());         // Pulsar Declination
  m_enphmin    = std::atof(parseParamList(params,3).c_str()); // minimum energy of extracted photons
  m_enphmax    = std::atof(parseParamList(params,4).c_str()); // minimum energy of extracted photons
  m_model      = std::atoi(parseParamList(params,5).c_str()); // choosen model
  m_seed       = std::atoi(parseParamList(params,6).c_str()); //Model Parameters: Random seed
  m_numpeaks   = std::atoi(parseParamList(params,7).c_str()); //Model Parameters: Number of peaks

 m_ppar1 = std::atof(parseParamList(params,8).c_str()); // model parameters
 m_ppar2 = std::atof(parseParamList(params,9).c_str());
 m_ppar3 = std::atof(parseParamList(params,10).c_str());
 m_ppar4 = std::atof(parseParamList(params,11).c_str());

  //Retrieve pulsar data from a list of DataList file.
  std::string pulsar_root = ::getenv("PULSARROOT");
  
  // if this is in the GLEAM environment, allow for separate path specified by env var PULSARDATA
  std::string pulsar_data(pulsar_root+"/data/"); // the default, perhaps overriden

  const char * gleam = ::getenv("PULSARDATA");
  if( gleam!=0) {
      pulsar_data =  std::string(gleam)+"/";
  }
  std::string ListFileName = pulsar_data + "PulsarDataList.txt";
  
  std::ifstream ListFile;
  ListFile.open(ListFileName.c_str(), std::ios::in);
  char DataListFileName[200];
  std::string CompletePathFileName = "";

  ListFile.getline(DataListFileName,200); 

  int Retrieved = 0;
  int AllRetrieved = 0;
 
  if (! ListFile.is_open()) 
    {
      std::cout << "Error opening file containing DataList files" << ListFileName  
		<< " (check whether $PULSARROOT is set" << std::endl; 
      exit(1);
    }  else 
      {

	ListFile >> DataListFileName;
	while (!ListFile.eof()) 
	  {	
	    CompletePathFileName = pulsar_data + std::string(DataListFileName);
	    Retrieved = getPulsarFromDataList(CompletePathFileName);
	    
	    if (Retrieved == 1)
	      {
		std::cout << "Pulsar " << m_PSRname << " found in file " << CompletePathFileName << std::endl;
		AllRetrieved = 1;
	      } 
	    ListFile >> DataListFileName ;	  

	  }

	if (AllRetrieved == 0)
	  {
	    std::cout << "Pulsar " << m_PSRname 
		      << " not found in any DataList file.Please Check if PULSARROOT is correctly set " 
		      << std::endl;
	    exit(1);
	  }
      }

  //Assign as starting ephemeris the first entry of the vectors... 
  m_t0Init = m_t0InitVect[0];
  m_t0 = m_t0Vect[0];
  m_t0End = m_t0EndVect[0];
  m_period = m_periodVect[0];
  m_pdot = m_pdotVect[0];
  m_p2dot = m_p2dotVect[0];
  m_f0 = m_f0Vect[0];
  m_f1 = m_f1Vect[0];
  m_f2 = m_f2Vect[0];
  m_phi0 = m_phi0Vect[0];


  //Init SolarSystem stuffs useful for barycentric decorrections

  astro::JulianDate JDStart(StartMissionDateMJD+JDminusMJD);
  m_earthOrbit = new astro::EarthOrbit(JDStart); 
  
  astro::SkyDir m_PulsarDir(m_RA,m_dec,astro::SkyDir::EQUATORIAL);
  m_PulsarVectDir = m_PulsarDir.dir();
  m_GalDir = std::make_pair(m_PulsarDir.l(),m_PulsarDir.b());
  m_l = m_GalDir.first;
  m_b = m_GalDir.second;


  //writes out an output log file

  std::string logLabel = m_PSRname + "Log.txt";
  ofstream PulsarLog(logLabel.c_str());
  
  PulsarLog << "\n********   PulsarSpectrum Log for pulsar" << m_PSRname << std::endl;
  PulsarLog << "**   Name : " << m_PSRname << std::endl;
  PulsarLog << "**   Position : (RA,Dec)=(" << m_RA << "," << m_dec 
	    << ") ; (l,b)=(" << m_l << "," << m_b << ")" << std::endl; 
  PulsarLog << "**   Flux above 100 MeV : " << m_flux << " ph/cm2/s " << std::endl;
  
  //Writes down on Log all the ephemerides
  for (unsigned int n=0; n < m_t0Vect.size(); n++)
    {
      PulsarLog << "**   Ephemerides valid from " << m_t0InitVect[n] 
		<< " to " << m_t0EndVect[n] << " (MJD): " << std::endl;
      PulsarLog << "**     Epoch (MJD) :  " << m_t0Vect[n] << std::endl;
      PulsarLog << std::setprecision(8) << "**     TxBary (MJD) where fiducial point (phi=0) is reached : " 
		<< m_txbaryVect[n] << std::endl;
      PulsarLog << "**     Phi0 (at Epoch t0) : " << m_phi0Vect[n] << std::endl;
      
      if (m_ephemType == "P")
	{
	  PulsarLog << "**     Ephemrides type: PERIOD" << std::endl;
	  PulsarLog << std::setprecision(14) << "**     Period : " << m_periodVect[n] << " s. | f0: " << m_f0Vect[n] << std::endl;
	  PulsarLog << std::setprecision(14) << "**     Pdot : " <<  m_pdotVect[n]  << " | f1: " << m_f1Vect[n] << std::endl; 
	  PulsarLog << std::setprecision(14) << "**     P2dot : " <<  m_p2dotVect[n]  << " | f2: " << m_f2Vect[n] << std::endl; 
	} 
      else if (m_ephemType == "F")
	{
	  PulsarLog << "**Ephemrides type: FREQUENCY" << std::endl;
	  PulsarLog << std::setprecision(14) << "**     Period : " << m_periodVect[n] << " s. | f0: " << m_f0Vect[n] << std::endl;
	  PulsarLog << std::setprecision(14) << "**     f1: " << m_f1Vect[n] << std::endl; 
	  PulsarLog << std::setprecision(14) << "**     f2: " << m_f2Vect[n] << std::endl; 
	}
    }

  PulsarLog << "**\n**   Enphmin : " << m_enphmin << " keV | Enphmax: " << m_enphmax << " keV" << std::endl;
  PulsarLog << "**   Mission started at (MJD) : " << StartMissionDateMJD << " (" 
	    << std::setprecision(12) << (StartMissionDateMJD+JDminusMJD)*SecsOneDay 
	    << " sec.)" << std::endl;
  PulsarLog << "**************************************************" << std::endl;

  //Instantiate an object of PulsarSim class
  m_Pulsar    = new PulsarSim(m_PSRname, m_seed, m_flux, m_enphmin, m_enphmax, m_period);
 
  //Instantiate an object of SpectObj class
  if (m_model == 1)
    {
      PulsarLog << "**   Model chosen : " << m_model << " --> Using Phenomenological Pulsar Model " << std::endl;  
      m_spectrum = new SpectObj(m_Pulsar->PSRPhenom(double(m_numpeaks), m_ppar1,m_ppar2,m_ppar3,m_ppar4),1);
      m_spectrum->SetAreaDetector(EventSource::totalArea());

      PulsarLog << "**   Effective Area set to : " << m_spectrum->GetAreaDetector() << " m2 " << std::endl; 

      if (DEBUG)
	{
	  std::cout << "**   Model chosen : " << m_model << " --> Using Phenomenological Pulsar Model " << std::endl;  
	  std::cout << "**  Effective Area set to : " << m_spectrum->GetAreaDetector() << " m2 " << std::endl; 
 	}  
    }
  else 
    {
      std::cout << "ERROR!  Model choice not implemented " << std::endl;
      exit(1);
    }


  PulsarLog.close();

  //Save the output txt file..
  int DbFlag = saveDbTxtFile();

  if (DEBUG) 
    if (DbFlag == 0)
      { 
	std::cout << "Database Output file created from scratch " << std::endl;
      } 
    else 
      {
	std::cout << "Appendended data to existing Database output file" << std::endl;
      }
}
/////////////////////////////////////////////////
PulsarSpectrum::~PulsarSpectrum() 
{  
  delete m_Pulsar;
  delete m_spectrum;
  delete m_earthOrbit;
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
 * This method find the time when the next photon will come. It takes into account decorrections due to
 * ephemerides and period derivatives and barycentric decorrections. This method also check for 
 * changes in the ephemerides validity range.
 * For the ephemerides decorrections the parameters used are:
 * <ul>
 * <li> <i>Epoch</i>, espressed in MJD;</li>
 * <li> <i>Initial Phase Phi0</i>;</li>
 * <li> <i>First Period derivative</i>;</li>
 * <li> <i>Second Period derivative</i>;</li>
 * </ul>
 * <br>
 * The barycentric decorrections takes into account conversion TDB->TT, Geometrical and Shapiro delays.
 * This method also calls the interval method in SpectObj class to determine the next photon in a rest-frame without 
 * ephemerides effects.
 */
double PulsarSpectrum::interval(double time)
{  

  double timeTildeDecorr = time + (StartMissionDateMJD)*SecsOneDay; //Arrival time decorrected
  double timeTilde = timeTildeDecorr + getBaryCorr(timeTildeDecorr); //should be corrected before applying ephem de-corrections
  
  double initTurns = getTurns(timeTilde); //Turns made at this time
  double intPart; //Integer part
  double tStart = modf(initTurns,&intPart)*m_period; // Start time for interval

  //Checks whether ephemerides (period,pdot, etc..) are within validity ranges
  if (((timeTilde/SecsOneDay) < m_t0Init) || ((timeTilde/SecsOneDay) > m_t0End)) 
    {
#if DEBUG
      std::cout << "Warning!Time is out of range of validity for pulsar " << m_PSRname 
		<< ": Switching to new ephemerides set..." << std::endl;
#endif
	for (int e=0; e < m_t0Vect.size();e++)
	  if (((timeTilde/SecsOneDay) > m_t0InitVect[e]) && ((timeTilde/SecsOneDay) < m_t0EndVect[e])) 
	    {
	      // std::cout << "\n" << timeTilde -(StartMissionDateMJD)*SecsOneDay <<  " PulsarSpectrum Phase is " 
	      // << std::setprecision(20) << getTurns(timeTilde) << "phi0 is " << m_phi0 << std::endl;

	      m_t0Init = m_t0InitVect[e];
	      m_t0 = m_t0Vect[e];
 	      m_t0End = m_t0EndVect[e];
	      m_f0 = m_f0Vect[e];
	      m_f1 = m_f1Vect[e];
	      m_f2 = m_f2Vect[e];
	      m_period = m_periodVect[e];
	      m_pdot = m_pdotVect[e];
	      m_p2dot = m_p2dotVect[e];
	      m_phi0 = m_phi0Vect[e];

	      if (DEBUG)
		{
		  std::cout << "Valid Ephemerides set found:" << std::endl;
		  std::cout << "MJD("<<m_t0Init<<"-"<<m_t0End<<") --> Epoch t0 = MJD " << m_t0 << std::endl; 
		  std::cout << "f0 " << m_f0 << " f1 " << m_f1 << " f2 " << m_f2 << std::endl;
		  std::cout << "Period " << m_period << " pdot " << m_pdot << " p2dot " << m_p2dot << std::endl;
		}

	      //Re-instantiate PulsarSim and SpectObj
	      delete m_Pulsar;
	      m_Pulsar = new PulsarSim(m_PSRname, m_seed, m_flux, m_enphmin, m_enphmax, m_period);

	      if (m_model == 1)
		{
		  delete m_spectrum;
		  m_spectrum = new SpectObj(m_Pulsar->PSRPhenom(double(m_numpeaks), m_ppar1,m_ppar2,m_ppar3,m_ppar4),1);
		  m_spectrum->SetAreaDetector(EventSource::totalArea());
		}

	      m_N0 = m_N0 + initTurns - getTurns(timeTilde); //Number of turns at next t0
	       double intN0;
	       double N0frac = modf(m_N0,&intN0); // Start time for interval
	       m_N0 = m_N0 - N0frac;
	       //	       std::cout << std::setprecision(20) << " Turns now are " << initTurns  
	       //              << " ; at t0 " << m_N0 << std::endl;	           
	       
	       //m_phi0 = m_phi0 - N0frac;
	       
	       if (DEBUG)
		 {
		   std::cout << std::setprecision(20) << " At Next t0 Number of Turns will be: " << m_N0 << std::endl;
		 }
	    }
	// else
	//	    {
	//std::cout << "Warning! Valid ephemerides not found!Proceeding with the old ones" << std::endl; 
	//}
    }

  if (DEBUG)
    {
      if ((int(timeTilde - (StartMissionDateMJD)*SecsOneDay) % 1000) < 1.5)
	std::cout << "**  Time reached is: " << timeTildeDecorr-(StartMissionDateMJD)*SecsOneDay
		  << " s. MET in TT, and " << timeTilde-(StartMissionDateMJD)*SecsOneDay
		  << " at the SSB (in TDB) ||| Pulsar " << m_PSRname << std::endl;
    }


   //Ephemerides calculations...

  double inte = m_spectrum->interval(tStart,m_enphmin); //deltaT (in system where Pdot = 0
  double inteTurns = inte/m_period; // conversion to # of turns
  double totTurns = initTurns + inteTurns + m_phi0; //Total turns at the nextTimeTilde
  double nextTimeTilde = retrieveNextTimeTilde(timeTilde, totTurns, (ephemCorrTol/m_period));


  double nextTimeDecorr = getDecorrectedTime(nextTimeTilde); //Barycentric decorrections
 
  if (DEBUG)
    { 
     std::cout << "\nTimeTildeDecorr at Spacecraft (TT) is: " 
	       << timeTildeDecorr - (StartMissionDateMJD)*SecsOneDay << "sec." << std::endl;
     std::cout << "TimeTilde at SSB (in TDB) is: " << timeTilde - (StartMissionDateMJD)*SecsOneDay << "sec." << std::endl;
     std::cout << "nextTimeTilde at SSB (in TDB) is:" << nextTimeTilde - (StartMissionDateMJD)*SecsOneDay << "sec." << std::endl;
     std::cout << " nextTimeTilde decorrected (TT)is" << nextTimeDecorr - (StartMissionDateMJD)*SecsOneDay << "sec." << std::endl;
     std::cout << " corrected is " << nextTimeDecorr + getBaryCorr(nextTimeDecorr) - (StartMissionDateMJD)*SecsOneDay << std::endl;
     std::cout << " interval is " <<  nextTimeDecorr - timeTildeDecorr << std::endl;
    }
  
  return nextTimeDecorr - timeTildeDecorr; //interval between the de-corected times
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
 * dN(t) = N<sub>0</sub> + phi<sub>0</sub> + f<sub>0</sub> + f<sub>1</sub>(t-t<sub>0</sub>) + 1/2f<sub>2</sub>(t-t<sub>0</sub>)<sup>2</sup>
 * </blockquote>
 * <br>where 
 * <blockquote>
 *<ul>
 * <li> t<sub>0</sub> is an ephemeris epoch.In this simulator epoch must be expressed in MJD;</li>
 * <li> N<sub>0</sub> is the number of turns at epoch t<sub>0</sub>;
 * <li> phi<sub>0</sub> is the phase shift at epoch t<sub>0</sub>;
 * <li>f<sub>0</sub> pulse frequency at time t<sub>0</sub> (usually given in Hz)</li>,
 * <li>f<sub>1</sub> the 1st time derivative of frequency at time t< sub>0</sub> (Hz s<sup>-1</sup>);</li>
 * <li>f<sub>2</sub> the 2nd time derivative of frequency at time t<sub>0</sub> (Hz s<sup>-2</sup>).</li>
 * </ul>
 * </blockquote>
 */

double PulsarSpectrum::getTurns( double time )
{
  double dt = time - m_t0*SecsOneDay;
  return m_phi0 + m_N0 + m_f0*dt + 0.5*m_f1*dt*dt + ((m_f2*dt*dt*dt)/6.0);
}

/////////////////////////////////////////////
/*!
 * \param tTilde Time from where retrieve the nextTime;
 * \param totalTurns Number of turns completed by the pulsar at nextTime;
 * \param err Phase tolerance between totalTurns and the number of turns at nextTime
 *
 * <br>In this method a recursive way is used to find the <i>nextTime</i>. Starting at <i>tTilde</i>, the method returns 
 * nextTime, the time where the number of turns completed by the pulsar is equal to totalTurns (within the choosen tolerance).  
 */
double PulsarSpectrum::retrieveNextTimeTilde( double tTilde, double totalTurns, double err )
{

  double tTildeUp = tTilde;
  double tTildeDown = tTilde;  
  double NTdown = totalTurns - getTurns(tTildeDown); 
  double NTup = totalTurns - getTurns(tTildeUp);
  int u = 0;
  while ((NTdown*NTup)>0)
    {
      u++;
      tTildeUp = tTilde+m_period*pow(2.0,u); 
      NTdown = totalTurns - getTurns(tTildeDown); 
      NTup = totalTurns - getTurns(tTildeUp);
    }
  
  double tTildeMid = (tTildeDown+tTildeUp)/2.0;
  double NTmid = totalTurns - getTurns(tTildeMid);

  while(fabs(NTmid) > err )
    { 
      if (fabs(NTmid) < err) break;
      NTmid = totalTurns - getTurns(tTildeMid);
      NTdown = totalTurns - getTurns(tTildeDown); 
      if ((NTmid*NTdown)>0)
	{
	  tTildeDown = tTildeMid;
	  tTildeMid = (tTildeDown+tTildeUp)/2.0;
	} 
      else 
	{
	  tTildeUp = tTildeMid;
	  tTildeMid = (tTildeDown+tTildeUp)/2.0;
	}
    }
  if (DEBUG)
    {
      std::cout << "**  RetrieveNextTimeTilde " << std::endl;
      std::cout << std::setprecision(30) << "  Stop up is " << tTildeUp << " NT " << NTup << std::endl;
      std::cout << std::setprecision(30) << "        down is " << tTildeDown << " NT " << NTdown << " u= " << u << std::endl;
      std::cout << std::setprecision(30) << "        mid is " << tTildeMid << " NT " << NTmid << " u= " << u << std::endl;
      std::cout << "     nextTimeTilde is " << tTildeMid << std::endl;
    }

  return tTildeMid;
}

/////////////////////////////////////////////
/*!
 * \param ttInput Photon arrival time at spacecraft in Terrestrial Time (expressed in MJD converted in seconds)
 *
 * <br>
 * This method computes the barycentric corrections for the photon arrival time at spacecraft and returns 
 * the time in Solar System Barycentric Time. The corrections implemented at the moment are:
 * <ul>
 *  <li> Conversion from TT to TDB;
 *  <li> Geometric Time Delay due to light propagation in Solar System;
 *  <li> Relativistic Shapiro delay;
 * </ul>   
 */
double PulsarSpectrum::getBaryCorr( double ttInput )
{
  //Start Date;
  astro::JulianDate ttJD(StartMissionDateMJD+JDminusMJD);
  ttJD = ttJD+(ttInput - (StartMissionDateMJD)*SecsOneDay)/SecsOneDay;

  if (DEBUG)
    {
      std::cout << std::setprecision(30) << "\nBarycentric Corrections for time " << ttJD << " (JD)" << std::endl;
    }

  // Conversion TT to TDB, from JDMissionStart (that glbary uses as MJDref)
  double tdb_min_tt = m_earthOrbit->tdb_minus_tt(ttJD);

  double timeMET = ttInput - (StartMissionDateMJD)*SecsOneDay;
  Hep3Vector scPos;

  //Exception error in case of time not in the range of available position (when using a FT2 file)

  try {
    scPos = astro::GPS::instance()->position(timeMET);
  }  catch (std::runtime_error  & eObj) {
    // check to see this is the exception I want to ignore, rethrowing if it is not:
    std::string message(eObj.what());
    if (message.find("WARNING: Time out of Range!") == std::string::npos) {
      throw;
    }
    // if so, do nothing....
  }

  //Correction due to geometric time delay of light propagation 
  Hep3Vector GeomVect = (scPos/clight) - m_solSys.getBarycenter(ttJD);
  double GeomCorr = GeomVect.dot(m_PulsarVectDir);

  //Correction due to Shapiro delay.
  Hep3Vector sunV = m_solSys.getSolarVector(ttJD);

  // Angle of source-sun-observer
  double costheta = - sunV.dot(m_PulsarVectDir) / ( sunV.mag() * m_PulsarVectDir.mag() );
  double m = 4.9271e-6; // m = G * Msun / c^3
  double ShapiroCorr = 2.0 * m * log(1+costheta);

  if (DEBUG)
    {
      std::cout << std::setprecision(20) << "** --> TDB-TT = " << tdb_min_tt << std::endl;
      std::cout << std::setprecision(20) << "** --> Geom. delay = " << GeomCorr << std::endl;
      std::cout << std::setprecision(20) << "** --> Shapiro delay = " << ShapiroCorr << std::endl;
      std::cout << std::setprecision(20) << "** ====> Total= " << tdb_min_tt + GeomCorr + ShapiroCorr  << " s." <<std::endl;
    }  

  return tdb_min_tt + GeomCorr + ShapiroCorr; //seconds

}

/////////////////////////////////////////////
/*!
 * \param CorrectedTime Photon arrival time at SSB (TDB expressed in MJD)
 *
 * <br>
 * This method returns the correspondent decorrected time starting from a photon arrival time
 * at the Solar System barycenter expressed in Barycentric Dynamical Time. This function uses the bisection method
 * for inverting barycentric corrections <br>
 * The corrections implemented at the moment are:
 * <ul>
 *  <li> Conversion from TT to TDB;
 *  <li> Geometric Time Delay due to light propagation in Solar System;
 *  <li> Relativistic Shapiro delay;
 * </ul>   
 */
double PulsarSpectrum::getDecorrectedTime(double CorrectedTime)
{
  //Use the bisection method to find the inverse of the de-corrected time, i.e. the time (in TT)
  //that after correction is equal to Time in TDB

  double deltaMax = 510.0; //max deltat (s)
  double ttUp = CorrectedTime + deltaMax;
  double ttDown = CorrectedTime - deltaMax;
  double ttMid = (ttUp + ttDown)/2;

  int i = 0;
  double hMid = 1e30; //for the 1st iteration

  while (fabs(hMid)>baryCorrTol )
    {
      double hDown = (CorrectedTime - (ttDown + getBaryCorr(ttDown)));
      hMid = (CorrectedTime - (ttMid + getBaryCorr(ttMid)));
                
      //std::cout << std::setprecision(30) 
      //	<< "\n" << i << "**ttUp " << ttUp << " ttDown " << ttDown << " mid " << ttMid << std::endl;
             
      if (fabs(hMid) < baryCorrTol) break;
      if ((hDown*hMid)>0)
	{
	  ttDown = ttMid;
	  ttMid = (ttUp+ttDown)/2;
	}
      else 
	{
	  ttUp = ttMid;
	  ttMid = (ttUp+ttDown)/2;
	}
       i++;
    }
   
   return ttMid;
}

/////////////////////////////////////////////
/*!
 * \param None
 *
 * <br>
 * This method gets from the ASCII <i>PulsarDatalist.txt</i> file stored in <i>/data</i> directory)
 * the names of the files that contain the pulsars parameters. The default one is <i>BasicDataList.txt</i>.
 * This method returns a integer status code (1 is Ok, 0 is failure)
 */
int PulsarSpectrum::getPulsarFromDataList(std::string sourceFileName)
{
  int Status = 0;
  std::ifstream PulsarDataTXT;
  
  if (DEBUG)
    {
      std::cout << "\nOpening Pulsar Datalist File : " << sourceFileName << std::endl;
    }
  
  PulsarDataTXT.open(sourceFileName.c_str(), std::ios::in);
  
  if (! PulsarDataTXT.is_open()) 
    {
      std::cout << "Error opening Datalist file " << sourceFileName  
		<< " (check whether $PULSARROOT is set" << std::endl; 
      Status = 0;
      exit(1);
    }
  
  char aLine[400];  
  PulsarDataTXT.getline(aLine,400); 

  char tempName[15] = "";

  double flux, ephem0, ephem1, ephem2, t0Init, t0, t0End, txbary, phi0, period, pdot, p2dot, f0, f1, f2, phas;
  char ephemType[15] = "";

  while (PulsarDataTXT.eof() != 1)
    {

      PulsarDataTXT >> tempName >> flux >> ephemType >> ephem0 >> ephem1 >> ephem2 >> t0Init >> t0 >> t0End  >> txbary;
     
      if (std::string(tempName) == m_PSRname)
	{

	  Status = 1;
	  m_flux = flux;
	  m_ephemType = ephemType;

	  m_t0InitVect.push_back(t0Init);
	  m_t0Vect.push_back(t0);
	  m_t0EndVect.push_back(t0End);
	  m_txbaryVect.push_back(txbary);

	  //Period-type ephemerides
	  if (std::string(ephemType) == "P")
	    {

	      period = ephem0;
	      pdot = ephem1;
	      p2dot = ephem2;
	      f0 = 1.0/period;
	      f1 = -pdot/(period*period);
	      f2 = 2*pow((pdot/period),2.0)/period - p2dot/(period*period);

	    } 
	  else if (std::string(ephemType) == "F")
	    {
	      //Frequency-style ephemrides
	      f0 = ephem0;
	      f1 = ephem1;
	      f2 = ephem2;
	      period = 1.0/f0;
	    }
	  
	  m_periodVect.push_back(period);
	  m_pdotVect.push_back(pdot);
	  m_p2dotVect.push_back(p2dot);
	  m_f0Vect.push_back(f0);
	  m_f1Vect.push_back(f1);
	  m_f2Vect.push_back(f2);

	  double dt = (txbary-t0)*SecsOneDay;

	  phi0 = -1.0*(f0*dt
		       + (f1/2.0)*dt*dt
		       + (f2/6.0)*dt*dt*dt); 
	
	  phi0 = modf(phi0,&phas);

	  if (phi0 < 0. ) 
	    phi0++;

	  m_phi0Vect.push_back(phi0);
	  
	}
    }

  return Status;
}

/////////////////////////////////////////////
/*!
 * \param None
 *
 * <br>
 * This method saves the relevant information in a file, named SimPulsar_spin.txt.
 * The format of this ASCII file is such that it can be given to gtpulsardb to produce
 * a D4-compatible FITS file, that can be used with pulsePhase
 */
int PulsarSpectrum::saveDbTxtFile()
{
  int Flag = 0;

  std::string DbOutputFileName = "SimPulsars_spin.txt";
  std::ifstream DbInputFile;

  //Checks if the file exists

  DbInputFile.open(DbOutputFileName.c_str(), std::ios::in);

  if (! DbInputFile.is_open()) 
    {
      Flag = 0;
    }
  else
    {
      Flag = 1;
    }

  DbInputFile.close();

  if (DEBUG)
    {
      std::cout << "Saving Pulsar ephemerides on file " << DbOutputFileName << std::endl;
    }

  std::ofstream DbOutputFile;

  if (Flag == 0)
    {
      DbOutputFile.open(DbOutputFileName.c_str(), std::ios::out);
      DbOutputFile << "# Simulated pulsars output file generated by PulsarSpectrum." << std::endl;
      DbOutputFile << "SPIN_PARAMETERS\n\n# Then, a column header."  << std::endl;
      DbOutputFile << "PSRNAME RA DEC EPOCH_INT EPOCH_FRAC TOAGEO_INT TOAGEO_FRAC TOABARY_INT TOABARY_FRAC ";
      DbOutputFile << "F0 F1 F2 RMS VALID_SINCE VALID_UNTIL BINARY_FLAG SOLAR_SYSTEM_EPHEMERIS OBSERVER_CODE" <<std::endl;
      DbOutputFile.close();
  }

  //Writes out the infos of the file
  DbOutputFile.open(DbOutputFileName.c_str(),std::ios::app);
  double tempInt, tempFract, f0, f1, f2;
  for (unsigned int ep = 0; ep < m_periodVect.size(); ep++)
    {
      DbOutputFile << "\"" << m_PSRname << std::setprecision(5) << "\" " << m_RA << " " << m_dec << " ";
      tempFract = modf(m_t0Vect[ep],&tempInt);
      DbOutputFile << std::setprecision (8) << tempInt << " " << tempFract << " ";
    
      tempFract = getDecorrectedTime(m_txbaryVect[ep]*SecsOneDay)/SecsOneDay;


      tempFract = modf(tempFract,&tempInt);
      DbOutputFile << std::setprecision (8) << tempInt << " " << tempFract << " ";

      tempFract = modf(m_txbaryVect[ep],&tempInt);
      DbOutputFile << std::setprecision (8) << tempInt << " " << tempFract << " ";

      DbOutputFile << std::setprecision(14) << m_f0Vect[ep] << " " << m_f1Vect[ep] << " " << m_f2Vect[ep] << " 0.0 " 
		   << m_t0InitVect[ep] << " " << m_t0EndVect[ep] 
      		   << " F \"JPL DE200\" P" << std::endl;

    }
  

  DbOutputFile.close();
  return Flag;
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
