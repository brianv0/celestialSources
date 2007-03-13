/////////////////////////////////////////////////
// File PulsarSpectrum.cxx
// Implementation of PulsarSpectrum class
//////////////////////////////////////////////////

#include "Pulsar/PulsarSpectrum.h"
#include "Pulsar/PulsarConstants.h"
#include "SpectObj/SpectObj.h"
#include "flux/SpectrumFactory.h"
#include <cmath>

#define DEBUG 0
#define BARYCORRLOG 0
#define BINDEMODLOG 0
#define TNOISELOG 0

using namespace cst;

//======================================================================
// copied from atErrors.h
//----------------------------------------------------------------------
/* Error Code of atFunctions. */
#define NOT_CONVERGED 10        /*equation was not solved*/
//----------------------------------------------------------------------
// copied from atKepler.c (& modified)
//----------------------------------------------------------------------
//#include "atFunctions.h"
//#include "atError.h"
//#include <math.h>

/*
 * solve Kepler equation (KEPLER)  g + e sin E = E
 */
int atKepler(
        double g,        /* input: mean anomaly */
        double eccent,        /* input: eccentricity */
        double *e)        /* output: eccentric anomaly */
{
    static double eps = 1e-7;
    static int imax = 50;

    int i;
    static double error, deltae, d__1;

    *e = g;
    if (g == 0.) return 0;

    for (i=0; i<imax; i++) {
        deltae = (g - *e + eccent * std::sin(*e)) / (1. - eccent * std::cos(*e));
        *e += deltae;
        error = (d__1 = deltae / *e, std::fabs(d__1));
        if (error < eps) return 0;
    }
    return NOT_CONVERGED;
}

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
  m_ff=false;
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
  m_f0NoNoise = 0.0;
  m_f1NoNoise = 0.0;
  m_f2NoNoise = 0.0;
  m_N0 = 0.0;
  m_t0Init = 0.0;
  m_t0 = 0.0;
  m_t0End = 0.0;  
  m_phi0 = 0.0;
  m_model = 0;
  m_seed = 0;
  m_TimingNoiseModel = 0;
  m_TimingNoiseActivity = 0.;
  m_TimingNoiseRMS = 0.;
  m_TimingNoiseMeanRate=0.;
  m_TimingNoiseTimeNextEvent=0.;
  m_BinaryFlag = 0;
  m_Porb = 0;
  m_Porb_dot = 0.;
  m_asini = 0;
  m_xdot = 0.;
  m_ecc = 0;
  m_ecc_dot = 0.;
  m_omega = 0;
  m_omega_dot = 0.;
  m_gamma = 0.;
  m_shapiro_r = 0.;
  m_shapiro_s = 0.;
  m_t0PeriastrMJD = 0;
  m_t0AscNodeMJD = 0;
  m_PPN =0.;

   //Read from XML file
  m_PSRname    = parseParamList(params,0).c_str();            // Pulsar name
  m_RA = std::atof(parseParamList(params,1).c_str());         // Pulsar Right Ascension
  m_dec = std::atof(parseParamList(params,2).c_str());         // Pulsar Declination
  m_enphmin    = std::atof(parseParamList(params,3).c_str()); // minimum energy of extracted photons
  m_enphmax    = std::atof(parseParamList(params,4).c_str()); // minimum energy of extracted photons
  m_model      = std::atoi(parseParamList(params,5).c_str()); // choosen model
  m_seed       = std::atoi(parseParamList(params,6).c_str()); //Model Parameters: Random seed

  //Setting random seed
  m_PSpectrumRandom = new TRandom;
  m_PSpectrumRandom->SetSeed(m_seed);
 
  if (m_model == 1) //Phenomenological model
    {
      m_ppar0   = std::atoi(parseParamList(params,7).c_str()); //Model Parameters: Number of peaks  
      m_ppar1 = std::atof(parseParamList(params,8).c_str());   // model parameters
      m_ppar2 = std::atof(parseParamList(params,9).c_str());
      m_ppar3 = std::atof(parseParamList(params,10).c_str());
      m_ppar4 = std::atof(parseParamList(params,11).c_str());
    }
  else 
    if (m_model == 2) //PulsarShape model
      {
	m_ppar0   = std::atoi(parseParamList(params,7).c_str());  //Model Parameters: Use normalization?
	m_PSRShapeName = parseParamList(params,8).c_str();        // model parameters
      }


  //Retrieve pulsar data from a list of DataList file.
  std::string pulsar_root = ::getenv("PULSARROOT");
  
  // if this is in the GLEAM environment, allow for separate path specified by env var PULSARDATA
  std::string pulsar_data(pulsar_root+"/data/"); // the default, perhaps overriden

  const char * gleam = ::getenv("PULSARDATA");
  if( gleam!=0) {
      pulsar_data =  std::string(gleam)+"/";
  }

  //Scan PulsarDataList.txt
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
		      << " not found in any DataList file.Please Check if $PULSARDATA is correctly set " 
		      << std::endl;
	    exit(1);
	  }
      }

  //  ListFile.getline(DataListFileName,200); 


  //Scan for Binary Pulsar if there are

  if (m_BinaryFlag ==1)
    {
      std::string BinListFileName = pulsar_data + "PulsarBinDataList.txt";
      
      //Look to BinDataList.txt for binary pulsar data
      std::ifstream BinListFile;
      BinListFile.open(BinListFileName.c_str(), std::ios::in);
      char BinDataListFileName[200];
      BinListFile.getline(BinDataListFileName,200); 

      if (! BinListFile.is_open()) 
	{
	  std::cout << "Error opening file containing Binary DataList files" << BinListFileName  
		    << " (check whether $PULSARDATA is set" << std::endl; 
	  exit(1);
	}  else 
	  {

	    BinListFile >> BinDataListFileName;
	    while (!BinListFile.eof()) 
	      {	
		CompletePathFileName = pulsar_data + std::string(BinDataListFileName);
		Retrieved = getOrbitalDataFromBinDataList(CompletePathFileName);
		
		if (Retrieved == 1)
		  {
		    std::cout << "Binary Pulsar " << m_PSRname << " found in file " << CompletePathFileName << std::endl;
		    AllRetrieved = 1;
		  } 
		BinListFile >> BinDataListFileName ;	  
	      }

	    if (AllRetrieved == 0)
	      {
		std::cout << "Binary Pulsar " << m_PSRname 
			  << " not found in any BinDataList file.Please Check if PULSARDATA is correctly set " 
			  << std::endl;
		exit(1);
	      }
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
  m_f0NoNoise = m_f0Vect[0];
  m_f1NoNoise = m_f1Vect[0];
  m_f2NoNoise = m_f2Vect[0];

  m_phi0 = m_phi0Vect[0];

  //Init SolarSystem stuffs useful for barycentric decorrections
  astro::JulianDate JDStart(StartMissionDateMJD+JDminusMJD);
  m_earthOrbit = new astro::EarthOrbit(JDStart); 
  
  astro::SkyDir m_PulsarDir(m_RA,m_dec,astro::SkyDir::EQUATORIAL);
  m_PulsarVectDir = m_PulsarDir.dir();
  m_GalDir = std::make_pair(m_PulsarDir.l(),m_PulsarDir.b());
  m_l = m_GalDir.first;
  m_b = m_GalDir.second;

  //Redirect output to a subdirectory
  const char * pulsarOutDir = ::getenv("PULSAROUTFILES");

  std::string logLabel;
  // override obssim if running in Gleam environment
  if( pulsarOutDir!=0) 
    logLabel = std::string(pulsarOutDir) + "/" + m_PSRname + "Log.txt";
  else
    logLabel = m_PSRname + "Log.txt";

  ofstream PulsarLog(logLabel.c_str());

  //Write infos to Log file  

  PulsarLog << "\n********   PulsarSpectrum Log for pulsar" << m_PSRname << std::endl;
  PulsarLog << "**   Name : " << m_PSRname << std::endl;
  PulsarLog << "**\n**   Position : (RA,Dec)=(" << m_RA << "," << m_dec 
	    << ") ; (l,b)=(" << m_l << "," << m_b << ")" << std::endl; 
  PulsarLog << "**\n**   Flux above 100 MeV : " << m_flux << " ph/cm2/s " << std::endl;
  PulsarLog << "**   Enphmin: " << m_enphmin << " keV | Enphmax: " << m_enphmax << " keV" << std::endl;
  PulsarLog << "**************************************************" << std::endl;
  
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
	  PulsarLog << "**     Ephemerides type: PERIOD" << std::endl;
	  PulsarLog << std::setprecision(14) << "**     Period : " << m_periodVect[n] << " s. | f0: " << m_f0Vect[n] << std::endl;
	  PulsarLog << std::setprecision(14) << "**     Pdot : " <<  m_pdotVect[n]  << " | f1: " << m_f1Vect[n] << std::endl; 
	  PulsarLog << std::setprecision(14) << "**     P2dot : " <<  m_p2dotVect[n]  << " | f2: " << m_f2Vect[n] << std::endl; 
	} 
      else if (m_ephemType == "F")
	{
	  PulsarLog << "**Ephemerides type: FREQUENCY" << std::endl;
	  PulsarLog << std::setprecision(14) << "**     Period : " << m_periodVect[n] << " s. | f0: " << m_f0Vect[n] << std::endl;
	  PulsarLog << std::setprecision(14) << "**     f1: " << m_f1Vect[n] << std::endl; 
	  PulsarLog << std::setprecision(14) << "**     f2: " << m_f2Vect[n] << std::endl; 
	}
    }

  PulsarLog << "**\n**   Mission Reference time: MJD " << StartMissionDateMJD << " (" 
	    << std::setprecision(12) << (StartMissionDateMJD+JDminusMJD)*SecsOneDay 
	    << " sec.)" << std::endl;

  if (m_TimingNoiseModel != 0)
    { 
      if (TNOISELOG==1)
	{
	  std::ofstream TimingNoiseLogFile((m_PSRname + "TimingNoise.log").c_str());
	  TimingNoiseLogFile << "Timing Noise log file using model: " << m_TimingNoiseModel << std::endl;
	  TimingNoiseLogFile << "tMET\tA\tf0_l\tf0_n\tf1_l\tf1_n\tf2_l\tf2_n\tPhiResid"<<std::endl;
	  TimingNoiseLogFile.close();
	}

      //define a default mean rate of about 1 day
      m_TimingNoiseMeanRate = 1/86400.;
      //interval according to Poisson statistics
      double startTime = Spectrum::startTime();
      m_TimingNoiseTimeNextEvent = startTime -log(1.-m_PSpectrumRandom->Uniform(1.0))/m_TimingNoiseMeanRate; 

      //Determine next Timing noise event according to the rate R m_TimingNoiseRate


      if (m_TimingNoiseModel == 1) // Timing model #1 - Delta8 parameter (Arzoumanian94)
	{
	  PulsarLog << "**\n**   Timing Noise Model : " << m_TimingNoiseModel << " (Stability parameter, Arzoumanian 1994)" << std::endl;
	  PulsarLog << "**      Timing Noise Events Mean Rate : " << m_TimingNoiseMeanRate << std::endl;
	}
      else if (m_TimingNoiseModel == 2) //Timing model #2 - of Random Walk (Cordes 1980) 
	{
	  PulsarLog << "**\n**   Timing Noise Model : " << m_TimingNoiseModel << " (Random Walk; Cordes 1980)" << std::endl;
	  PulsarLog << "**      Timing Noise Events Mean Rate : " << m_TimingNoiseMeanRate << std::endl;
	  double alpha = 0.57;
	  double PdotCrab = 422E-15;
	  double SigmaRMSCrab = 0.012/0.033;
	  m_TimingNoiseActivity = -1.0*alpha*std::log10(PdotCrab)+alpha*std::log10(m_pdotVect[0]);
	  m_TimingNoiseRMS = pow(10.,m_TimingNoiseActivity)*SigmaRMSCrab;
	  PulsarLog << "**          Cordes Activity Parameter A: " <<  m_TimingNoiseActivity << std::endl;
	  PulsarLog << "**          Cordes Activity Timing Noise RMS (phase): " <<  m_TimingNoiseRMS << std::endl;
	}
    }

  PulsarLog << "**************************************************" << std::endl;

  //Add lines in case of binary

  if (m_BinaryFlag == 1)
    {
      PulsarLog << "**   Pulsar in a Binary System! Orbital Data:" << std::endl;
      PulsarLog << "**     Orbital period: " << m_Porb << " s." << std::endl;
      PulsarLog << "**     Projected major semiaxis (a * sini): " << m_asini << " lightsec." <<std::endl; 
      PulsarLog << "**     Eccentricity: " << m_ecc << std::endl;
      PulsarLog << "**     Longitude of periastron: " <<  m_omega << " deg." << std::endl;
      PulsarLog << "**     Epoch of Periastron (MJD): " << m_t0PeriastrMJD << std::endl;
      PulsarLog << "**     Epoch of Ascending Node (MJD): " << m_t0AscNodeMJD << std::endl;
      if (m_PPN ==0)
	PulsarLog << "**   No Post Newtonian Parameterization " << std::endl; 

      PulsarLog << "**************************************************" << std::endl;

      
      if ((m_BinaryFlag ==1) && (BINDEMODLOG))
 	{
	  std::ofstream BinDemodLogFile((m_PSRname + "BinDemod.log").c_str());
	  BinDemodLogFile << "tMET\tdt\tE\tAt\tOmega\tecc\tasini\tdt_Roemer\tdt_einstein\tdt_shapiro" << std::endl;
	  BinDemodLogFile.close();
	}
   
      if (BARYCORRLOG)
	{
	  std::ofstream BaryCorrLogFile((m_PSRname + "BaryCorr.log").c_str());
	  BaryCorrLogFile << "\nPulsar" << m_PSRname << " : Barycentric corrections log generated by PulsarSpectrum" << std::endl;
	  BaryCorrLogFile << "tMET\tTDB-TT\tGeomDelay\tShapiroDelay" << std::endl;
	  BaryCorrLogFile.close();
	}
    }

  //Instantiate an object of PulsarSim class
  m_Pulsar    = new PulsarSim(m_PSRname, m_seed, m_flux, m_enphmin, m_enphmax, m_period);
 
  //Instantiate an object of SpectObj class
  if (m_model == 1)
    {
      PulsarLog << "**   Model chosen : " << m_model << " --> Using Phenomenological Pulsar Model " << std::endl;  
      m_spectrum = new SpectObj(m_Pulsar->PSRPhenom(double(m_ppar0), m_ppar1,m_ppar2,m_ppar3,m_ppar4),1);
      m_spectrum->SetAreaDetector(EventSource::totalArea());

      PulsarLog << "**   Effective Area set to : " << m_spectrum->GetAreaDetector() << " m2 " << std::endl; 

      if (DEBUG)
	{
	  std::cout << "**   Model chosen : " << m_model << " --> Using Phenomenological Pulsar Model " << std::endl;  
	  std::cout << "**  Effective Area set to : " << m_spectrum->GetAreaDetector() << " m2 " << std::endl; 
 	}  
    }
  else 
    if (m_model == 2)
      {
	PulsarLog << "**   Model chosen : " << m_model << " --> Using External 2-D Pulsar Shape" << std::endl;  
	m_spectrum = new SpectObj(m_Pulsar->PSRShape(m_PSRShapeName,m_ppar0),1);
	m_spectrum->SetAreaDetector(EventSource::totalArea());
	
	PulsarLog << "**   Effective Area set to : " << m_spectrum->GetAreaDetector() << " m2 " << std::endl; 
	
	if (DEBUG)
	  {
	    std::cout << "**   Model chosen : " << m_model << " --> Using External 2-D Pulsar Shape" << std::endl;  
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


  //Save the binary data output txt file..
  if (m_BinaryFlag ==1)
    {

      int BinDbFlag = saveBinDbTxtFile();
      
      if (DEBUG) 
	if (BinDbFlag == 0)
	  { 
	    std::cout << "Database for Binary pulsars file created from scratch " << std::endl;
	  } 
	else 
	  {
	    std::cout << "Database for Binary pulsars appended to existing binary Database output file" << std::endl;
	  }
    }

}

/////////////////////////////////////////////////
PulsarSpectrum::~PulsarSpectrum() 
{  
  delete m_Pulsar;
  delete m_spectrum;
  delete m_earthOrbit;
  delete m_PSpectrumRandom;
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
  //this should be corrected before applying barycentryc decorr + ephem de-corrections
  double timeTilde = timeTildeDecorr + getBaryCorr(timeTildeDecorr,0); 

  double timeTildeDemodulated =0.;
  
  //binary demodulation
  if (m_BinaryFlag ==0)
    {
      timeTildeDemodulated = timeTilde;
    }
  else
    {
      timeTildeDemodulated = getIterativeDemodulatedTime(timeTilde,1);      
    }

  //Phase assignment



  double intPart=0.; //Integer part
  double PhaseNoNoise,PhaseWithNoise=0.;

  //pippo

  //Apply timing noise
  if (m_TimingNoiseModel !=0)
    {
      //Check for a next Timing noise event
      if ((timeTildeDemodulated-(StartMissionDateMJD)*SecsOneDay) > m_TimingNoiseTimeNextEvent)
	{
	  //interval according to Poisson statistics
	  m_TimingNoiseTimeNextEvent+= -log(1.-m_PSpectrumRandom->Uniform(1.0))/m_TimingNoiseMeanRate;
 
	  if (DEBUG)
	    {
	      std::cout << std::setprecision(30)<<"Timing Noise Event!Next Event at t=" 
			<< m_TimingNoiseTimeNextEvent << " |dt=" 
			<< m_TimingNoiseTimeNextEvent
		-(timeTildeDemodulated-(StartMissionDateMJD)*SecsOneDay) <<std::endl;
	    }

	  if (m_TimingNoiseModel ==1) // Timing Noise #1
	    {
	      m_f2 = 0;
	      PhaseNoNoise = getTurns(timeTildeDemodulated);
	      PhaseNoNoise = modf(PhaseNoNoise,&intPart); // phase for linear evolution 
	      if ( PhaseNoNoise <0.)
		PhaseNoNoise+=1.;
	      m_TimingNoiseActivity = 6.6 + 0.6*log10(m_pdot) + m_PSpectrumRandom->Gaus(0,0.5);
	      
	      //estimate an f2
	      double Sign = m_PSpectrumRandom->Uniform();
	      if (Sign > 0.5)
		m_f2 = ((m_f0*6.*std::pow(10.,m_TimingNoiseActivity))*1e-24);
	      else
		m_f2 = -1.0*((m_f0*6.*std::pow(10.,m_TimingNoiseActivity))*1e-24);
	    }
	  else  if (m_TimingNoiseModel == 2) // Timing Noise #2
	    {
	      m_f2 = 0.;
	      double tempPhi0 = m_phi0; 
	      m_phi0 = 0.;
	      PhaseNoNoise = getTurns(timeTildeDemodulated);
	      PhaseNoNoise = modf(PhaseNoNoise,&intPart); // phase for linear evolution 
	      if ( PhaseNoNoise <0.)
		PhaseNoNoise+=1.;

	      m_phi0 = tempPhi0+m_PSpectrumRandom->Gaus(0,m_TimingNoiseRMS);  
	    }

	  
	  PhaseWithNoise = getTurns(timeTildeDemodulated);
	  PhaseWithNoise = modf(PhaseWithNoise,&intPart); // phase for linear evolution 
	  if ( PhaseWithNoise <0.)
	    PhaseWithNoise+=1.;

	  
	  std::ofstream TimingNoiseLogFile((m_PSRname + "TimingNoise.log").c_str(),std::ios::app);
	  m_f2NoNoise = 0.;
	  double ft_l = GetFt(timeTildeDemodulated,m_f0NoNoise,m_f1NoNoise,m_f2NoNoise);
	  double ft_n = GetFt(timeTildeDemodulated,m_f0,m_f1,m_f2);
	  double ft1_l = GetF1t(timeTildeDemodulated,m_f1NoNoise,m_f2NoNoise);
	  double ft1_n = GetF1t(timeTildeDemodulated,m_f1,m_f2);
	  double ft2_l = m_f2NoNoise;//
	  double ft2_n = m_f2;//
	 
	  if (TNOISELOG==1)
	    {
	      TimingNoiseLogFile << std::setprecision(30) << timeTildeDemodulated-(StartMissionDateMJD)*SecsOneDay
				 << "\t" << m_TimingNoiseActivity
				 << "\t" <<ft_l << "\t" << ft_n 
				 << "\t" <<ft1_l << "\t" << ft1_n 
				 << "\t" <<ft2_l << "\t" << ft2_n 			 
				 << "\t" << PhaseNoNoise << "\t" << PhaseWithNoise << std::endl;
	    }
	  

	  if (DEBUG)
	    {
	      std::cout << std::setprecision(30) << " Activity=" << m_TimingNoiseActivity 
			<< "f2=" << m_f2 << " PN=" << PhaseWithNoise 
			<< " dPhi=" << PhaseWithNoise-PhaseNoNoise <<std::endl;
	    }
	}
    }

  /*
  if (m_TimingNoiseModel ==1)
    {
	  m_f2 = 0.;
	  initTurnsNoNoise = getTurns(timeTildeDemodulated);
	  
	  Delta8 = 6.6 + 0.6*log10(m_pdot) + m_PSpectrumRandom->Gaus(0,1.0);
	  
	  double Sign = m_PSpectrumRandom->Uniform();
	  if (Sign > 0.5)
	    m_f2 = ((m_f0*6.*std::pow(10.,Delta8))*1e-24);
	  else
	    m_f2 = -1.0*((m_f0*6.*std::pow(10.,Delta8))*1e-24);
    }
  else if (m_TimingNoiseModel ==2)
    {

      if ((timeTildeDemodulated-(StartMissionDateMJD)*SecsOneDay) > m_timeMETNextNoiseEvent)
	{
	  m_timeMETNextNoiseEvent+=400.;
	  double sigma_r_Crab_1628 = 0.012; // rms residual for crab
	  double sigma_r = std::sqrt(sigma_r_Crab_1628*pow(10.,m_TimingNoiseActivity)); 
	  double sigma_m = 8e-5;
	  double sigma_R = sigma_r;
	  double coeff = 15.5;
	  // m_RWStrength_FN = coeff*coeff*(((sigma_r*sigma_r)-(sigma_m*sigma_m))/(sigma_R*sigma_R));
	  std::cout << std::setprecision(30) << (timeTildeDemodulated-(StartMissionDateMJD)*SecsOneDay)
		    << "Timing Noise Event! next will be at " << m_timeMETNextNoiseEvent << std::endl;
	  //" s. MET rms:" << sigma_r 
	  //    << " Strenght: " << m_RWStrength_FN << std::endl;
	
	  int CurrentDay = int((timeTildeDemodulated-(StartMissionDateMJD)*SecsOneDay)/86400);
	  double CurrentDayFrac = 86400*(((timeTildeDemodulated-(StartMissionDateMJD)*SecsOneDay)/86400)-CurrentDay);
	  //std::cout << std::setprecision(30) << CurrentDay << "+ " << CurrentDayFrac << std::endl;
	  
	  //  if (CurrentDayFrac < 100.)
	  //{
	  m_RWRate = 1/86400.;
	  //std::cout << "Current day " << CurrentDay << std::endl;
	  double DeltaPhi0 = m_PSpectrumRandom->Gaus(0.,std::sqrt((m_RWStrength_PN/m_RWRate)));
	  double Deltaf0 = m_PSpectrumRandom->Gaus(0.,std::sqrt((m_RWStrength_FN/m_RWRate)));
	  double Deltaf1 = m_PSpectrumRandom->Gaus(0.,std::sqrt((m_RWStrength_SN/m_RWRate)));
	  std::cout << std::setprecision(30) << "dPhi0=" << DeltaPhi0 << " df0=" << Deltaf0 
		    << " df1=" << Deltaf1 << " f2" << m_f2 << std::endl;

	  m_phi0 = m_phi0 +  DeltaPhi0;
	  m_f0 = m_f0 + Deltaf0;
	  m_f1 = m_f1 + Deltaf1;
	  m_f2 = 0.;
	  
	  //std::cout << std::setprecision(30) << "Current day " << CurrentDayFrac 
	  //	    << " phi0:" << m_phi0 << " f0:" << m_f0 << " nonoise" << m_f0NoNoise << " df0" << Deltaf0 << std::endl;
	  
	  //}
	 }
    }
  */


  double initTurns = getTurns(timeTildeDemodulated); //Turns made at this time
  double tStart = modf(initTurns,&intPart)*m_period; // Start time for interval
  if (tStart < 0.)
    tStart+=m_period;

  //  double DeltaPhiNoise = modf(initTurns,&intPart)-modf(initTurnsNoNoise,&intPart);

  if (DEBUG)
    {
      std::cout << std::setprecision(30) << "\n" << timeTilde -(StartMissionDateMJD)*SecsOneDay 
		<<  " turns are " << getTurns(timeTilde) 
		<<  " phase is " << tStart/m_period << " phi0 is " << m_phi0 << std::endl;
    }

  //Checks whether ephemerides (period,pdot, etc..) are within validity ranges
  if (((timeTildeDemodulated/SecsOneDay) < m_t0Init) || ((timeTildeDemodulated/SecsOneDay) > m_t0End)) 
    {
#if DEBUG
      std::cout << "Warning!Time is out of range of validity for pulsar " << m_PSRname 
		<< ": Switching to new ephemerides set..." << std::endl;
#endif
	for (unsigned int e=0; e < m_t0Vect.size();e++)
	  if (((timeTildeDemodulated/SecsOneDay) > m_t0InitVect[e]) && ((timeTildeDemodulated/SecsOneDay) < m_t0EndVect[e])) 
	    {
	
	      m_t0Init = m_t0InitVect[e];
	      m_t0 = m_t0Vect[e];
 	      m_t0End = m_t0EndVect[e];
	      m_f0 = m_f0Vect[e];
	      m_f1 = m_f1Vect[e];
	      m_f2 = m_f2Vect[e];
	      m_f0NoNoise = m_f0Vect[e];
	      m_f1NoNoise = m_f1Vect[e];
	      m_f2NoNoise = m_f2Vect[e];
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
		  m_spectrum = new SpectObj(m_Pulsar->PSRPhenom(double(m_ppar0), m_ppar1,m_ppar2,m_ppar3,m_ppar4),1);
		  m_spectrum->SetAreaDetector(EventSource::totalArea());
	      
		}

	      m_N0 = m_N0 + initTurns - getTurns(timeTildeDemodulated); //Number of turns at next t0
	      double intN0;
	      double N0frac = modf(m_N0,&intN0); // Start time for interval
	      m_N0 = m_N0 - N0frac;
	      std::cout << std::setprecision(20) << " Turns now are " << initTurns  
			<< " ; at t0 " << m_N0 << std::endl;	           
	       
	      //m_phi0 = m_phi0 - N0frac;
	       
	      if (DEBUG)
		{
		  std::cout << std::setprecision(20) << " At Next t0 Number of Turns will be: " << m_N0 << std::endl;
		}
	    }
	 else
	   {
	     std::cout << "Warning! Valid ephemerides not found!Proceeding with the old ones" << std::endl; 
	   }
    }

  if (DEBUG)
    {
      if ((int(timeTilde - (StartMissionDateMJD)*SecsOneDay) % 1000) < 1.5)
	std::cout << "\n\n** " << m_PSRname << " Time is: " 
		  << timeTildeDecorr-(StartMissionDateMJD)*SecsOneDay << " s MET in TT | "
		  << timeTilde-(StartMissionDateMJD)*SecsOneDay << "s. MET in SSB in TDB (barycorr.) | "
		  << timeTildeDemodulated-(StartMissionDateMJD)*SecsOneDay << "s. MET in SSB in TDB (barycorr. + demod.)| "
		  << std::endl;
    }

  //Log out the barycentric corrections  
  double bl=0;
  if (BARYCORRLOG)
     bl = getBaryCorr(timeTilde,1);
  
  //Ephemerides calculations...

  double inte = m_spectrum->interval(tStart,m_enphmin); //deltaT (in system where Pdot = 0
  double inteTurns = inte/m_period; // conversion to # of turns
  double totTurns = initTurns + inteTurns + m_phi0; //Total turns at the nextTimeTilde

  //Applying timingnoise	     
  
  double nextTimeTildeDemodulated = retrieveNextTimeTilde(timeTildeDemodulated, totTurns, (ephemCorrTol/m_period));
  
  /*
  if ((m_TimingNoiseModel != 0) && (TNOISELOG==1))
    {
      std::ofstream TimingNoiseLogFile((m_PSRname + "TimingNoise.log").c_str(),std::ios::app);
      m_f2NoNoise = 0.;
      double ft_l = GetFt(timeTildeDemodulated,m_f0NoNoise,m_f1NoNoise,m_f2NoNoise);
      double ft_n = GetFt(timeTildeDemodulated,m_f0,m_f1,m_f2);
      double ft1_l = 12.;//
      double ft1_n = 23;//
      double ft2_l = 11.;//
      double ft2_n = 33;//
      

      TimingNoiseLogFile << std::setprecision(30) << nextTimeTildeDemodulated-(StartMissionDateMJD)*SecsOneDay
			 << "\t" << m_TimingNoiseActivity
			 << "\t" <<ft_l << "\t" << ft_n 
			 << "\t" <<ft1_l << "\t" << ft1_n 
			 << "\t" <<ft2_l << "\t" << ft2_n 			 
			 << "\t" << (tStart/m_period)-PhaseNoise << std::endl;


   }
  */

  double nextTimeTilde = 0.;
  double nextTimeTildeDecorr = 0.;

  //inverse of binary demodulation and of barycentric corrections
  if (m_BinaryFlag == 0)
    {
      nextTimeTilde = nextTimeTildeDemodulated;
      nextTimeTildeDecorr = getDecorrectedTime(nextTimeTilde); //Barycentric decorrections
    }
  else 
    {
      nextTimeTilde = getBinaryDemodulationInverse(nextTimeTildeDemodulated);
      nextTimeTildeDecorr = getDecorrectedTime(nextTimeTilde); //Barycentric decorrections
    }
 
  if (DEBUG)
    { 
      std::cout << std::setprecision(15) << "\nTimeTildeDecorr at Spacecraft (TT) is: " 
		<< timeTildeDecorr - (StartMissionDateMJD)*SecsOneDay << "sec." << std::endl;
      std::cout << "Arrival time after start of the simulation : " 
		<< timeTildeDecorr - (StartMissionDateMJD)*SecsOneDay -Spectrum::startTime() << " s." << std::endl; 
      std::cout << std::setprecision(15) <<"  TimeTilde at SSB (in TDB) is: " 
		<< timeTilde - (StartMissionDateMJD)*SecsOneDay << "sec." << std::endl;
      std::cout << std::setprecision(15) <<"  TimeTildeDemodulated : " 
		<< timeTildeDemodulated - (StartMissionDateMJD)*SecsOneDay 
		<< "sec.; phase=" << modf(initTurns,&intPart) << std::endl;
      std::cout << std::setprecision(15) <<"  nextTimeTildeDemodulated at SSB (in TDB) is:" 
		<< nextTimeTildeDemodulated - (StartMissionDateMJD)*SecsOneDay << "s." << std::endl;
      std::cout << std::setprecision(15) <<"  nextTimeTilde at SSB (in TDB) is:" 
		<< nextTimeTilde - (StartMissionDateMJD)*SecsOneDay << "s." << std::endl;
      std::cout << std::setprecision(15) <<"  nextTimeTilde decorrected (TT)is" 
		<< nextTimeTildeDecorr - (StartMissionDateMJD)*SecsOneDay << "s." << std::endl;
      std::cout << " corrected is " 
		<< nextTimeTildeDecorr + getBaryCorr(nextTimeTildeDecorr,0) - (StartMissionDateMJD)*SecsOneDay << std::endl;
      std::cout << std::setprecision(15) <<"  -->Interval is " <<  nextTimeTildeDecorr - timeTildeDecorr << std::endl;
    }
  
  double interv = nextTimeTildeDecorr - timeTildeDecorr;
  if (interv < 0.)
    interv = m_periodVect[0]/100.;
  
  return interv;
  
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

  double tTildeUp,NTup = 0.;
  double tTildeDown,NTdown = 0;  

  tTildeUp = tTilde;
  tTildeDown = tTilde;  

  NTdown = totalTurns - getTurns(tTildeDown); 
  NTup = totalTurns - getTurns(tTildeUp);


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
double PulsarSpectrum::getBaryCorr( double ttInput, int LogCorrFlag)
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
  CLHEP::Hep3Vector scPos;

  //Exception error in case of time not in the range of available position (when using a FT2 file)

  try {
    //    std::cout<<astro::GPS::time()<<std::endl;
    //astro::GPS::update(timeMET);
    //    std::cout<<astro::GPS::time()<<std::endl;
    astro::GPS::instance()->time(timeMET);
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
  //GLAST position
  CLHEP::Hep3Vector GeomVect = (scPos/clight);
  double GLASTPosGeom = GeomVect.dot(m_PulsarVectDir);

  //Earth position
  GeomVect = - m_solSys.getBarycenter(ttJD);
  double EarthPosGeom = GeomVect.dot(m_PulsarVectDir);

  GeomVect = (scPos/clight) - m_solSys.getBarycenter(ttJD);

  double GeomCorr = GeomVect.dot(m_PulsarVectDir);

  //Correction due to Shapiro delay.
  CLHEP::Hep3Vector sunV = m_solSys.getSolarVector(ttJD);

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


  if ((LogCorrFlag == 1) && (ttInput-StartMissionDateMJD*SecsOneDay > 0.))
    {
	  std::ofstream BaryCorrLogFile((m_PSRname + "BaryCorr.log").c_str(),std::ios::app);
	  BaryCorrLogFile << std::setprecision(20)
			  << timeMET << "\t" <<  GLASTPosGeom << "\t"<< EarthPosGeom << "\t" 
			  << GeomCorr << "\t" << tdb_min_tt << "\t" << ShapiroCorr << std::endl;
	  BaryCorrLogFile.close();
    }

  return tdb_min_tt + GeomCorr + ShapiroCorr; //seconds
  
}


/////////////////////////////////////////////
/*!
 * \param tInput Photon arrival time do be demodulated
 * \param LogFlag Flag for writing output
 *
 * <br>
 * This method compute the binary demodulation in an iterative way 
 * using the getBinaryDemodulation method
 */
double PulsarSpectrum::getIterativeDemodulatedTime(double tInput, int LogFlag)
{

  double timeDemodulated = tInput+getBinaryDemodulation(tInput,0);
            
  double TargetTime = tInput;
  double delay = getBinaryDemodulation(tInput,0);

  int i=0;
  while ((fabs(timeDemodulated-delay-TargetTime) > DemodTol) && ( i < 50))
    {

      i++;
      timeDemodulated = TargetTime+delay;
      delay = getBinaryDemodulation(timeDemodulated,0);
    }
  
  if (LogFlag == 1)
    {
      delay = getBinaryDemodulation(timeDemodulated,1);
    }
  
  return timeDemodulated;
}

/////////////////////////////////////////////
/*!
 * \param tInput Photon arrival time do be demodulated
 * \param LogFlag Flag for writing output
 *
 * <br>
 * This method compute the binary demodulation using the orbital parameters
 * of the pulsar. The corrections computed are:
 * 
 *<ul>
 * <li> Roemer delay
 * <li> Einstein delay
 * <li> Shapiro delay
 *</ul>
 */
double PulsarSpectrum::getBinaryDemodulation( double tInput, int LogDemodFlag)
{
  double BinaryRoemerDelay = 0.;
  double dt = tInput-m_t0PeriastrMJD*SecsOneDay;
  double OmegaMean = 2*M_PI/m_Porb;
  double EccAnConst = OmegaMean*(dt - 0.5*(dt*dt*(m_Porb_dot/m_Porb)));
  
  //Calculate Eccenctric Anomaly solving Kepler equation using the atKepler function
  double EccentricAnomaly = 0.; 
  int status = atKepler(EccAnConst, m_ecc, &EccentricAnomaly); //AtKepler
  
  // atKepler not converged
  if (0 != status) {
     throw std::runtime_error("atKepler did not converge.");
  } 
  
  //Calculate True Anomaly
  double TrueAnomaly = 2.0 * std::atan(std::sqrt((1.0+m_ecc)/(1.0-m_ecc))*std::tan(EccentricAnomaly*0.5));
  TrueAnomaly = TrueAnomaly + 2*M_PI*floor((EccentricAnomaly - TrueAnomaly)/ (2*M_PI));
  while ((TrueAnomaly - EccentricAnomaly) > M_PI) TrueAnomaly -= 2*M_PI;
  while ((EccentricAnomaly - TrueAnomaly) > M_PI) TrueAnomaly += 2*M_PI;
  
  // Calculate longitude of periastron using m_omega_dot
  double Omega = DegToRad*m_omega + DegToRad*m_omega_dot*(TrueAnomaly/OmegaMean);
  
  // compute projected semimajor axis using x_dot
  double asini = m_asini + m_xdot*dt;
  
  //Binary Roemer delay
  BinaryRoemerDelay = ((std::cos(EccentricAnomaly)-m_ecc)*std::sin(Omega)
			     + std::sin(EccentricAnomaly)*std::cos(Omega)*std::sqrt(1-m_ecc*m_ecc));
  //Einstein Delay
  double BinaryEinsteinDelay = m_gamma*std::sin(EccentricAnomaly);
 
  //Shapiro binary delay
  double BinaryShapiroDelay = -2.0*m_shapiro_r*log(1.-m_ecc*std::cos(EccentricAnomaly)-m_shapiro_s*BinaryRoemerDelay);

  BinaryRoemerDelay = asini*BinaryRoemerDelay;
     
  if (DEBUG)
    {
      std::cout << "\n**  Binary modulations t=" << std::setprecision(15) << tInput-(StartMissionDateMJD)*SecsOneDay
		<< " dt=" << tInput - m_t0PeriastrMJD*SecsOneDay << std::endl;
      std::cout << "**  Ecc.Anomaly=" << EccentricAnomaly << " deg. True Anomaly=" << TrueAnomaly << " deg." << std::endl;
      std::cout << "**  Omega=" << Omega << " rad. Major Semiaxis " << asini << " light-sec" << std::endl;
      std::cout << "**  --> Binary Roemer Delay=" << BinaryRoemerDelay << " s." << std::endl;
      std::cout << "**  --> Binary Einstein Delay=" << BinaryEinsteinDelay << " s." << std::endl;
      std::cout << "**  --> Binary Shapiro Delay=" << BinaryShapiroDelay << " s." << std::endl;
  }
  
  if ((LogDemodFlag ==1) && (tInput-StartMissionDateMJD*SecsOneDay > 0.))
    {
      std::ofstream BinDemodLogFile((m_PSRname + "BinDemod.log").c_str(),std::ios::app);
      BinDemodLogFile << std::setprecision(15) 
		      << tInput-StartMissionDateMJD*SecsOneDay << "\t" << tInput - m_t0PeriastrMJD*SecsOneDay << "\t"
		      << EccentricAnomaly << "\t" << TrueAnomaly << "\t"
		      << Omega << "\t" << m_ecc << "\t" << asini << "\t"
		      << BinaryRoemerDelay << "\t" 
		      << BinaryEinsteinDelay << "\t" 
		      << BinaryShapiroDelay << "\t" << std::endl;
      BinDemodLogFile.close();
    }

  return -1.*(BinaryRoemerDelay+BinaryEinsteinDelay+BinaryShapiroDelay);
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
      double hDown = (CorrectedTime - (ttDown + getBaryCorr(ttDown,0)));
      hMid = (CorrectedTime - (ttMid + getBaryCorr(ttMid,0)));
                
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
 * \param CorrectedTime Photon arrival time Modulatedat SSB (TDB expressed in MJD)
 *
 * <br>
 * This method returns the correspondent inverse-demodulated time of each photons, includingRoemer delay, Einstein delay 
 * and Shapiro delay
 */
double PulsarSpectrum::getBinaryDemodulationInverse( double CorrectedTime)
{

  //Set initial conditions
  double deltaMax = m_asini*(1-m_ecc)+m_asini*std::sqrt(1-m_ecc*m_ecc); //max deltat (s)
  double deltaMin = m_asini*(1+m_ecc)-m_asini*std::sqrt(1-m_ecc*m_ecc);

  double tModUp = CorrectedTime + deltaMax;
  double tModDown = CorrectedTime - deltaMin;
  double tModMid = 0.;

  int i = 0;
  double hModUp = 0.;
  double hModDown = 0.;
  double hModMid = 1e30; //for the 1st iteration

   while (fabs(hModMid)>InverseDemodTol )
    {
      i++;
      tModMid = tModUp*0.5 + tModDown*0.5;
      
      if (tModMid == tModDown)
	{
	  tModDown = tModDown-1e-5*m_PSpectrumRandom->Uniform();
	  if (DEBUG) std::cout << "\n**Correcting down to tDown:" << tModDown << std::endl;
	  tModMid = tModUp*0.5 + tModDown*0.5;
	}

      if (tModMid == tModUp)
	{
	  tModUp = tModUp+1e-5*m_PSpectrumRandom->Uniform();
	  if (DEBUG) std::cout << "\n**Correcting up to tUp:" << tModUp << std::endl;
	  tModMid = tModUp*0.5 + tModDown*0.5;
	}

      hModMid = CorrectedTime - (getIterativeDemodulatedTime(tModMid,0));

      /*
       std::cout << std::setprecision(30) << "\n****Rand" << RandFrac << " delta " << (tModUp - tModDown)*RandFrac << std::endl;
       std::cout << std::setprecision(30) << " tUp-tmid " << tModUp - tModMid << std::endl;
       std::cout << std::setprecision(30) << " tDown-tmid " << tModDown - tModMid << std::endl;
       std::cout << std::setprecision(30) << "   tdown=" << tModDown << " tup=" << tModUp << " tMid=" << tModMid << std::endl;
       std::cout << std::setprecision(30) << "   hdown=" << hModDown << " hup=" << hModUp << " hMid=" << hModMid << std::endl;
      */
      
      if ((hModDown*hModMid)<0)
	{
	  tModUp = tModMid;
	  //std::cout << " tMid--->tTup ::: tUp=" << tModUp<<std::endl; 
	}
      else
	{
	  tModDown = tModMid;
	  //	  std::cout << " tMid--->tTDown ::: tDown=" << tModDown<<std::endl; 
	}

      hModUp = CorrectedTime - (getIterativeDemodulatedTime(tModUp,0));
      hModDown = CorrectedTime - (getIterativeDemodulatedTime(tModDown,0));

      if (i==50) 
	{
	  if (DEBUG)
	    {
	      std::cout << " ERROR!!!: Inverse demodulation does not converge! " << std::endl;
	    }
	  break;
    	}
    }
   
  return tModMid;

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

  double startTime = Spectrum::startTime();
  double endTime = astro::GPS::instance()->endTime();

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
  int tnmodel, binflag;

  while (PulsarDataTXT.eof() != 1)
    {

      PulsarDataTXT >> tempName >> flux >> ephemType >> ephem0 >> ephem1 >> ephem2 
		    >> t0Init >> t0 >> t0End  >> txbary >> tnmodel >> binflag;
      
      if (std::string(tempName) == m_PSRname)
	{

	  Status = 1;
	  m_flux = flux;
	  m_TimingNoiseModel = tnmodel;
	  m_ephemType = ephemType;
	  m_BinaryFlag = binflag;

	  //Check if txbary or t0 are before start of the simulation
	  double startMJD = StartMissionDateMJD+(startTime/86400.)+(550/86400.);
	  // std::cout << "T0 " << t0 << " start " << startMJD << std::endl;
	  if (t0 < startMJD)
	    {
	      std::cout << std::setprecision(10) << "Warning! Epoch t0 out the simulation range (t0-tStart=" << (t0-startMJD)
			<< " s.: changing to MJD " << startMJD << std::endl;
	      t0 = startMJD;
	    }
	  
	  if (txbary < startMJD)
	    {
	      std::cout << "Warning! txbary out the simulation range (t0-tStart=" << (txbary-startMJD)
			<< " s.: changing to MJD " << startMJD << std::endl;
	      txbary = startMJD;
	    }

	  //Check if txbary or t0 are after start of the simulation
	  double endMJD = StartMissionDateMJD+(endTime/86400.)-(550/86400.);
	  // std::cout << "T0 " << t0 << " start " << startMJD << std::endl;
	  if (t0 > endMJD)
	    {
	      std::cout << std::setprecision(10) << "Warning! Epoch t0 out the simulation range (t0-tEnd=" << (t0-endMJD)
			<< " s.: changing to MJD " << endMJD << std::endl;
	      std::cout << "***end at " << endTime << " corresp to " << endMJD << std::endl;
	      t0 = endMJD;
	    }
	  
	  if (txbary > endMJD)
	    {
	      std::cout << "Warning! txbary out the simulation range (t0-tEnd=" << (txbary-endMJD)
			<< " s.: changing to MJD " << endMJD << std::endl;
	      txbary = endMJD;
	    }

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
 * This method gets the orbital parameters of the binary pulsar from a 
 * <i>BinDataList/i> file among those specified in the<i>/data/PulsarBinDataList</i> 
 *
 * The orbital parameters are:
 *
 * m_Porb Orbital period;
 * m_asini Projected major semiaxis;
 * m_eccentr eccentricity;
 * m_omega longitude of periastron;
 * m_t0PeriastronMJD  epoch of periastron passage;
 * m_t0AscNodeMJD Epoch of the ascending node;
 *
 * The key for finding pulsar is the name of the files that contain the pulsars parameters. 
 * The default one is <i>BasicDataList.txt</i>.
 * Extra parameters are used to specify a PPN parameterization for General Relativity;
 * This method returns a integer status code (1 is Ok, 0 is failure)
 */
int PulsarSpectrum::getOrbitalDataFromBinDataList(std::string sourceBinFileName)
{
  int Status = 0;
  std::ifstream PulsarBinDataTXT;
  
  if (DEBUG)
    {
      std::cout << "\nOpening Binary Pulsar BinDatalist File : " << sourceBinFileName << std::endl;
    }
  
  PulsarBinDataTXT.open(sourceBinFileName.c_str(), std::ios::in);
  
  if (! PulsarBinDataTXT.is_open()) 
    {
      std::cout << "Error opening BinDatalist file " << sourceBinFileName  
		<< " (check whether $PULSARDATA is set)" << std::endl; 
      Status = 0;
      exit(1);
    }
  
  char aLine[400];  
  PulsarBinDataTXT.getline(aLine,400); 
 
  char tempName[30];
  double porb,asini,ecc,omega,t0peri,t0asc,ppn;

  while (PulsarBinDataTXT.eof() != 1)
    {
      
      PulsarBinDataTXT >> tempName >> porb >> asini >> ecc >> omega >> t0peri >> t0asc >> ppn;
      
      if (std::string(tempName) == m_PSRname)
	{
	  m_Porb = porb;
	  m_asini = asini;
	  m_ecc = ecc;
	  m_omega = omega;
	  m_t0PeriastrMJD = t0peri;
	  m_t0AscNodeMJD = t0asc;
	  m_PPN = ppn;
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
  double tempInt, tempFract;
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
/*!
 * \param None
 *
 * <br>
 * This method saves the relevant orbital parameters in a file, named SimPulsar_bin.txt.
 * The format of this ASCII file is such that it can be given to gtpulsardb to produce
 * a D4-compatible FITS file, that can be used with pulsePhase
 */
int PulsarSpectrum::saveBinDbTxtFile()
{
  int Flag = 0;

  std::string DbBinOutputFileName = "SimPulsars_bin.txt";
  std::ifstream DbBinInputFile;

  //Checks if the file exists

  DbBinInputFile.open(DbBinOutputFileName.c_str(), std::ios::in);

  if (! DbBinInputFile.is_open()) 
    {
      Flag = 0;
    }
  else
    {
      Flag = 1;
    }

  DbBinInputFile.close();

  if (DEBUG)
    {
      std::cout << "Saving Pulsar Orbital Data on file " << DbBinOutputFileName << std::endl;
    }

  std::ofstream DbBinOutputFile;

  if (Flag == 0)
    {
      DbBinOutputFile.open(DbBinOutputFileName.c_str(), std::ios::out);
      DbBinOutputFile << "# Simulated pulsars orbital data output file generated by PulsarSpectrum." << std::endl;
      DbBinOutputFile << "ORBITAL_PARAMETERS\n\n# This file can be converted to a D4 fits file using:"  << std::endl;
      DbBinOutputFile << "# >gtpulsardb SimPulsars_bin.txt" << std::endl;
      DbBinOutputFile << "PSRNAME PB PBDOT A1 XDOT ECC ECCDOT OM OMDOT T0 GAMMA OBSERVER_CODE SOLAR_SYSTEM_EPHEMERIS SHAPIRO_R SHAPIRO_S" << std::endl;;
      DbBinOutputFile.close();
  }

  //Writes out the infos of the file
  DbBinOutputFile.open(DbBinOutputFileName.c_str(),std::ios::app);
  
  DbBinOutputFile << "\"" << m_PSRname << "\" ";                                   // pulsar name
  DbBinOutputFile << std::setprecision(10) << m_Porb << " " << m_Porb_dot << " ";   // Orbital period and derivative
  DbBinOutputFile << std::setprecision(10) << m_asini << " " << m_xdot << " ";      // Projected semi-mayor axis and derivative
  DbBinOutputFile << std::setprecision(10) << m_ecc << " " << m_ecc_dot << " ";     // Eccentricity and derivative
  DbBinOutputFile << std::setprecision(10) << m_omega << " " << m_omega_dot << " "; // Long. Periastron and derivative
  DbBinOutputFile << std::setprecision(10) << m_t0PeriastrMJD << " " <<  m_gamma << " "; // t0 of Periastron and PPN gamma   
  DbBinOutputFile << "MR \"JPL DE200\" ";                                          // Observer code and ephemerides
  DbBinOutputFile << std::setprecision(10) << m_shapiro_r << " " << m_shapiro_s << std::endl; //Shapiro Parameters

  DbBinOutputFile.close();


  //In this case a summary D4 file is created
  std::ofstream DbSumInputFile("SimPulsars_summary.txt");
  DbSumInputFile << "SimPulsars_spin.txt\nSimPulsars_bin.txt" <<std::endl;
  DbSumInputFile.close();

  return Flag;

}


//////////////////////////////////////////////////////////
// no longer used
/////////////////////////////////////////////////////////
double PulsarSpectrum::GetEccentricAnomaly(double mytime)
{
  double OmegaMean = 2*M_PI/m_Porb;
  double dtime = (mytime-m_t0PeriastrMJD*SecsOneDay);
  double EccAnConst = OmegaMean*(dtime - 0.5*(m_Porb_dot/m_Porb)*dtime*dtime);

  double Edown = 0.;
  double EccAnDown = Edown-(m_ecc*std::sin(Edown))-EccAnConst;

  double Eup = 2*M_PI;
  double EccAnUp = Eup-(m_ecc*std::sin(Eup))-EccAnConst;

  double Emid = 0.5*(Eup + Edown);
  double EccAnMid = Emid-(m_ecc*std::sin(Emid))-EccAnConst;

  int i=0;
  while (fabs(EccAnMid) > 5e-7)
    {
      if ((EccAnDown*EccAnMid) < 0)
	{
	  Eup = Emid;
	  EccAnUp = EccAnMid;
	}
      else
	{
	  Edown = Emid;
	  EccAnDown = EccAnMid;
	}
      i++;
      Emid = 0.5*(Eup + Edown);
      EccAnMid = Emid-(m_ecc*std::sin(Emid))-EccAnConst;
      if (fabs(EccAnMid) < 5e-7)
	break;
    }

  return Emid;
}


/////////////////////////////////////////////
double PulsarSpectrum::energy(double time)
{
  return m_spectrum->energy(time,m_enphmin)*1.0e-3; //MeV
}

/////////////////////////////////////////////
/*!
 * \param time input time for calculating f(t)
 *
 * This method computes the frequency at a given instant t
 *
 */ 
double PulsarSpectrum::GetFt(double time, double myf0, double myf1, double myf2)
{
  double dt = time - m_t0*SecsOneDay;
  return myf0 + myf1*dt + 0.5*myf2*dt*dt;
}

/////////////////////////////////////////////
/*!
 * \param time input time for calculating f'(t)
 *
 * This method computes the frequency first derivativeat a given instant t
 *
 */ 
double PulsarSpectrum::GetF1t(double time, double myf1, double myf2)
{
  double dt = time - m_t0*SecsOneDay;
  return myf1 + myf2*dt;
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
