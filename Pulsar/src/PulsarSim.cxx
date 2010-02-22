//////////////////////////////////////////////////
// File PulsarSim.cxx
// contains the code for the implementation of the models
//////////////////////////////////////////////////

#include <cstdlib>
#include <iostream>
#include "Pulsar/PulsarConstants.h"
#include "Pulsar/PulsarSim.h"
#include "facilities/commonUtilities.h"

#define DEBUG 0
#define SAVETIMEPROFILE 0

using namespace cst;

//////////////////////////////////////////////////
/*!
 * \param name Name of the pulsar; 
 * \param seed Seed of the random generator;
 * \param flux Flux of the pulsar in ph/cm2/s; 
 * \param enphmin Minimun energy of the extracted photons in keV;
 * \param enphmax Maximum energy of the extracted photons in keV;
 * \param period Period of the pulsar;
*/
PulsarSim::PulsarSim(std::string name, int seed, double flux, double enphmin, double enphmax, double period)
{

  m_name = name;         //Pulsar name
  m_seed = seed;         //Random seed
  m_flux = flux;         //ph/cm2/s
  m_enphmin = enphmin;   // KeV
  m_enphmax = enphmax;   //KeV
  m_period  = period;    //s
  m_Tbin = Tbin;         //If a template is specified, m_Tbin depends upon the the bins contained in it.

  //look for pulsar data directory, following, in the order $PULSARDATA, $SKYMODEL_DIR/pulsars or $PULSARROOT/data
  //  std::string m_pulsardata_dir;
  
  if (::getenv("PULSARDATA"))
    {
      const char * psrdata = ::getenv("PULSARDATA");
      m_pulsardata_dir = std::string(psrdata);
    }
  else if (::getenv("SKYMODEL_DIR"))
    {
      const char * psrdata = ::getenv("SKYMODEL_DIR");
      m_pulsardata_dir =  std::string(psrdata)+"/pulsars";
    }
  else
    {
      const char * psrdata = ::getenv("PULSARROOT");
      m_pulsardata_dir =  std::string(psrdata)+"/data";
    }

  if (DEBUG)
    {
      std::cout << "PULSARDATA used is: "<< m_pulsardata_dir << std::endl;
    }


  if (::getenv("PULSAR_OUTPUT_LEVEL"))
    {
      const char * outlevel = ::getenv("PULSAR_OUTPUT_LEVEL");
      m_OutputLevel = atoi(outlevel);
      if ((m_OutputLevel<0) || (m_OutputLevel>2))
	m_OutputLevel=1;
    }
  else
    {
      m_OutputLevel=1;
    }


  if (m_OutputLevel > 0)
    {
      //Redirect output to a subdirectory
      const char * pulsarOutDir = ::getenv("PULSAROUTFILES");
      
      // override obssim if running in Gleam environment
      if( pulsarOutDir!=0) 
	m_LogFileName = std::string(pulsarOutDir) + "/" + m_name + "Log.txt";
      else
	m_LogFileName = m_name + "Log.txt";
      
      char temp[200];
      sprintf(temp,"**  Output Level set to: %d",m_OutputLevel);
      WriteToLog(std::string(temp));
    }




}

//////////////////////////////////////////////////
/*!
 * \param par0 Number of peaks(1) : 1 - Only one peak;
 *                                  2 - Two peaks;
 *                                  3 - LightCurve from an user-defined Time Profile txt file (1)
 * \param par1 Parameter En expressed in GeV;
 * \param par2 Parameter E0 expressed in GeV;
 * \param par3 Parameter g ;
 * \param par4 Parameter b ;
 *
 *(1) - <b>Note</b>: If you want to use the option 3 to obtain the lightcurve from a template, you should 
 *      have the txt file in the <i>data</i> directory.For Example if you have in the PulsarDataList a pulsar
 *      named PSRTEST you must have a correspondant file PSRTESTTimeProfile.txt in the <i>data</I. directory.<br>
 * This method creates a ROOT TH2D histogram according to a phenomenological model based on observations of known
 * gamma-ray pulsars. The 2D histogram is obtained by multiplying the lightcurve and the spectrum. 
 * The lightcurve is obtained by random generating a profile of 1 or 2 peaks separated by a minimum distance,
 * (that currently is set to one half of the period).
 * Otherwise you can have the lightcurve starting from a TXT file that contain the time profile. In this case
 * you should have the txt file in the <i>data</i> directory.For Example if you have in the PulsarDataList a pulsar
 * named PSRTEST you must have a correspondant file PSRTESTTimeProfile.txt in the <i>data</i>. directory.
 * <br>The spectrum is generated from an analytical form :<br>
 *
 * \image html NJSnForm.jpg
 * <br>
 * where the parameters E0,En,g and g are the parameters of the method. This form describes a power law spectrum
 * with an exponential cutoff. This formula is choosen according to Nel, De Jager (1995,see Ref.below ) and 
 * De Jager (2002, see Ref.below). According to Cheng (1994, see Ref. below) an Outer Gap scenario can be obtained 
 * by setting b=1.
 * The constant K indicates the normalisation. We choose to set the normalisation in accord to the 3rd EGRET Catalog
 * where the fluxes are reported above 100MeV in ph/s/cm2. 
 * The ROOT histogram is then saved to a ROOT file with the same name of the pulsar, and also a Txt Time profile is
 * saved. <br>
 * <i>For more informations and for a brief tutorial please see</i>:
 * <br>
 * <a href="#dejager02">http://www.pi.infn.it/~razzano/Pulsar/PulsarSpTutor/PulsarSpTutor.htm</a>
 * <br><br>
 * <i><b>References:</b></i>
 * <ul>
 *  <li><a name="chengpaper94"></a>Cheng, K.S. and Ding, W.K.Y.:1994, <i>Astrophys.J.</i>431,724;</li>
 *  <li><a name="neilpaper95"></a>Neil, H.I. and De Jager, O.C.:1995, <i>Astrophysics and Space Science</i> 230:209-306;</li>
 *  <li><a name="dejager02"></a>De Jager, O.C.:2002, <i>African Skies</i>,No 7;</li>
 *  <li><a name="egret3cat"></a>Hartman, R.C. et al.: 1999, <i>The Astrophysical Journal Supplement Series</i>,123:79-202;</li>
 * </ul>
 */
TH2D* PulsarSim::PSRPhenom(double par0, double par1, double par2, double par3, double par4)
{

  m_numpeaks = int(par0);

  // Part 1 - Spectrum

  double En = par1;
  double G1 = par3;
  double E0 = par2;
  double b =  par4;
  double K1 = 138e-8; //This is the constant, will be overwritten by normalization

  //Establish the lower and upper limits for the energy in ROOT histogram.Write out these infos
  double LowEnBound = TMath::Min(cst::EnNormMin,m_enphmin); 
  if (m_enphmax < cst::EnNormMax)
    m_enphmax = cst::EnNormMax;
  double HighEnBound = TMath::Max(cst::EnNormMax,m_enphmax); 

  //Writes out informations about the model parameters on the log file
  WriteToLog("******** Pulsar Phenomenological Model ********");

  char temp[300];
  sprintf(temp,"**  Random seed for the model : %d",m_seed);
  WriteToLog(std::string(temp));
  WriteToLog("**  Spectrum parameters: ");
  sprintf(temp,"**           En: %.2e | E0: %.2E",En,E0);
  WriteToLog(std::string(temp));
  sprintf(temp,"**           G1: %.2f | b: %.2f",G1,b);  
  WriteToLog(std::string(temp));
  sprintf(temp,"**  enphmin: %.2e keV | enphmax: %.2e keV",m_enphmin,m_enphmax);
  WriteToLog(std::string(temp));
  sprintf(temp,"**           Normalisation between %.2e keV and %.2e keV ",cst::EnNormMin,cst::EnNormMax);
  WriteToLog(std::string(temp));
  sprintf(temp,"**           Photon extraction between %.2e keV and %.2e keV ",m_enphmin,m_enphmax);
  WriteToLog(std::string(temp));
  sprintf(temp,"**  Spectrum calculated between %.2e keV and %.2e keV",LowEnBound,HighEnBound);
  WriteToLog(std::string(temp));

  //Create the spectrum profile
  double de = pow(HighEnBound/LowEnBound,1.0/Ebin);
 
  TF1 PulsarSpectralShape("PulsarSpectralShape", 
			  "([0]*((x/[1])^[2])*exp(-1.0*((x/[3])^[4])))", LowEnBound, HighEnBound);
  PulsarSpectralShape.SetParameters(K1,En,G1,E0,b);

  TF1 PulsarTimeCurve("PulsarTimeCurve",
		      "([2]*(1/(((x-[0])^2)+(([1]/2)^2))) + [5]*(1/(((x-[3])^2)+(([4]/2)^2))))",
		      0, m_period);

  TH1D TimeProfileLightCurve;//("TimeProfileLightCurve","TimeProfileLightCurve",m_Tbin,0,m_period);

  // Part 2- LightCurve

  double dt = 0.0;

  if ((m_numpeaks == 1) || (m_numpeaks == 2)) //case of random Lorentz peak generation
    {
      dt = m_period/(m_Tbin-1);

      double fwhm1,fwhm2,peak1,peak2,ampl1,ampl2,mindist=0.;;
      mindist = 0.5*m_period; //minimum distance in seconds (minPhase * m_period)
      
      //Set the Random engine.
      TRandom engine;
      engine.SetSeed(m_seed);
      
      // LightCurve generation:
      ampl1 = engine.Uniform();
      while (ampl1 < 0.1)  ampl1 = engine.Uniform();
      
      ampl2 = engine.Uniform();
      while (ampl2 < 0.1)  ampl2 = engine.Uniform();
      
      peak1 = engine.Uniform()*m_period; 
      fwhm1 = engine.Uniform()*m_period;
      
      while ((peak1 > 0.5*m_period)
	     || (peak1 < (fwhm1*2)))
	{
	  peak1 = engine.Uniform()*m_period; 
	  fwhm1 = engine.Uniform()*m_period;
	}    
      
      peak2 = engine.Uniform()*m_period; 
      fwhm2 = engine.Uniform()*m_period;
      
      while ((peak2 > (m_period-2*fwhm2))
	     || (peak2 <(peak1+mindist)))
	{
	  peak2 = engine.Uniform()*m_period; 
	  fwhm2 = engine.Uniform()*m_period;
	}

      if (fwhm1 < (0.01*m_period))
	fwhm1 = 0.01*m_period;
      
      if (fwhm2 < (0.01*m_period))
	fwhm2 = 0.01*m_period;
      
      
      //Remove first or second peak.
      if (m_numpeaks == 1)
	{
	  if (engine.Uniform() < 0.5) 
	    {
	      ampl1 = 0;
	    }
	  else 
	    {
	      ampl2 = 0;
	    }
	}

      char temp[200];
      printf(temp,"**\n**  Lightcurve parameters: (mindist = %.2f s.)",mindist);
      WriteToLog(std::string(temp));
      
      if (ampl1 !=0)
	{
	  sprintf(temp,"**           Peak 1 t = %.3f (ph= %.3f), FWHM: %.3f, ampl: %.3f",
		  peak1,peak1/m_period,fwhm1,ampl1);
	  WriteToLog(std::string(temp));
	  WriteToLog("***********************************************");
	}
      if (ampl2 !=0)
	{
	  sprintf(temp,"**           Peak 1 t = %.3f (ph= %.3f), FWHM: %.3f, ampl: %.3f",
		  peak2,peak2/m_period,fwhm2,ampl2);
	  WriteToLog(std::string(temp));
	  WriteToLog("***********************************************");
	}
      

      //    PulsarTimeCurve.SetParameters(peak1,fwhm1,ampl1,peak2,fwhm2,ampl2);  
  
    }
  else if (m_numpeaks == 3)
    {

      //Look for NAMETimeProfile.txt in the data directory...

      std::string TimeProfileFileName = facilities::commonUtilities::joinPath(m_pulsardata_dir, m_name+"TimeProfile.txt");

      // override obssim if running in Gleam environment
      //     if( gleam!=0) TimeProfileFileName = facilities::commonUtilities::joinPath(std::string(gleam), m_name + "TimeProfile.txt");

      if (m_OutputLevel>1)
	{
	  WriteToLog("Building lightcurve from file "+ TimeProfileFileName);
	}
      
      ifstream TimeProfileFile;
      TimeProfileFile.open(TimeProfileFileName.c_str(), std::ios::in);
  
      if (! TimeProfileFile.is_open()) 
	{
	  std::cout << "Error opening TimeProfile file " << TimeProfileFileName << " for Pulsar " << m_name << std::endl;
	  exit (1);
	}

      std::vector <double> timeCounts;
      double tempCount = 0.0;
      int b = 0;
      while (TimeProfileFile.eof() != 1)
	{
	  TimeProfileFile >> b >> tempCount;
	  timeCounts.push_back(tempCount);
	}

      TH1D TempProfile("TempProfile","TempProfile",timeCounts.size()-1,0,m_period);
      
      for (unsigned int i =0; i < timeCounts.size()-1; i++)
	{
	  TempProfile.SetBinContent(i,timeCounts[i]);
	}
      TimeProfileLightCurve = TempProfile;
      m_Tbin = TimeProfileLightCurve.GetNbinsX();
      dt = m_period/(m_Tbin-1);
    }

  // Part 3 - Combination of Spectrum and lightCurve and filling of TH2D

  double *e = new double[Ebin +1];
  for(int i = 0; i<=Ebin; i++)
    {
      e[i] = LowEnBound*pow(de,1.0*i); //KeV
    }

  gDirectory->Delete("Nv");

  m_Nv = new TH2D("Nv","Nv",m_Tbin,0.,m_period,Ebin, e);

  //Filling the TH2D Histogram
  double t = 0.0;
  for(int ti = 0; ti<m_Tbin; ti++)
    {
      t = ti*dt;
      double nt = 0.0;
      if ((m_numpeaks == 1) || (m_numpeaks == 2))
	{
	  nt = PulsarTimeCurve.Eval(t);
	}
      else if (m_numpeaks == 3 )
	{
	  nt = TimeProfileLightCurve.GetBinContent(ti);
	}

      for(int ei = 0; ei < Ebin; ei++)
	{
	  double nv = PulsarSpectralShape.Eval(e[ei]);
	  m_Nv->SetBinContent(ti+1, ei+1, nt*nv);// [ph/(cm² s KeV)]
	}
    }

  // !!!  Conversion 1/cm² -> 1/m² IMPORTANT m_Nv has to be in [ph/m2s KeV)
  // !!!  BUT in the XML file the flux is espressed in ph/cm2/s according to EGRET Catalogs

  m_Nv->Scale(1.0e+4);  // [ph/(m² s keV)]
  m_flux*= 1.0e+4;      //  ph/m²/s
  
  // nph = nv * dE * dt
  
  TH2D *nph = Nph(m_Nv); //ph/m²
  
  int ei2 = nph->GetYaxis()->FindBin(cst::EnNormMin);
  int ei3 = nph->GetYaxis()->FindBin(cst::EnNormMax);
 
  //Normalisation factor according to band between EGRET1 and EGRET2 energies
  //Integration is on a averaged flux over period
  double norm = m_Nv->Integral(0,m_Tbin,ei2,ei3,"width")/m_period; // ph/m2/s

  m_Nv->Scale(m_flux/norm);

  delete nph;
  delete[] e;
  
  if (::getenv("PULSAR_OUT_NV"))
    {
      SaveNv(m_Nv); // ph/m2/s/keV
    }

  if (SAVETIMEPROFILE)
    {
      SaveTimeProfile(m_Nv); //ph/m2/s/kev
    }

  return m_Nv;
  delete m_Nv;
}


//////////////////////////////////////////////////
/*!
 * \param ModelShapeName Name of the file containing model
 * \param NormalizeFlux 0 - Use the normalization of the TH2D histogram
 *                      1 - Normalize the TH2D histogram to the m_flux value (in ph/cm2/s above 100 MeV)
 * 
 * This method use an external TH2D histogram containing an arbitrary phase-energy distribution and 
 * use it to extract photons, instead of using an analytical form for the spectrum as in PSRPhenom method
 * The model named by the parameter ModelShapeName must be located in the PULSARDATA directory
 */
TH2D* PulsarSim::PSRShape(std::string ModelShapeName, int NormalizeFlux)
{
  //  std::cout << "Testing new model Shape for pulsar " << m_name << " using normalization option " << NormalizeFlux << std::endl; 

  //Writes out an output log file for the pulsar
 
  //Redirect output to a subdirectory
  const char * pulsarOutDir = ::getenv("PULSAROUTFILES");

  std::string logSimLabel;
  // override obssim if running in Gleam environment
  if( pulsarOutDir!=0) 
    logSimLabel = facilities::commonUtilities::joinPath(std::string(pulsarOutDir), m_name + "Log.txt");
  else
    logSimLabel = m_name + "Log.txt"; 

  //TimeProfileFileName = std::string(gleam)+"/"+ m_name + "TimeProfile.txt";

  //Writes out informations about the model parameters on the file
  WriteToLog("******** Pulsar Shape Model ********");
  WriteToLog("** Using and arbitrary 2-d shape as spectrum");
  WriteToLog("** Pulsar Shape used: "+ std::string(ModelShapeName));
  if (NormalizeFlux ==1)
    WriteToLog("** Use normalization :");
  else
    WriteToLog("** No normalization used:");
 

  //Look for ModelShapeName.root in the data directory...
  std::string ModelShapeInputFileName = facilities::commonUtilities::joinPath(m_pulsardata_dir, ModelShapeName + ".root");

  WriteToLog("** Using Shape "+std::string(ModelShapeName)+ " located at: "+std::string(ModelShapeInputFileName));

  //Load 2HD Histogram

  TFile *ModelShapeInputFile = new TFile(ModelShapeInputFileName.c_str());
  TH2D *m_Nv = (TH2D*) ModelShapeInputFile->Get("Nv"); //ph m^(-2) s^(-1) keV^(-1)
    
  //  double shapePhmin = m_Nv->GetXaxis()->GetXmin();
  //double shapePhmax = m_Nv->GetXaxis()->GetXmax();
  double shapeEmin = m_Nv->GetYaxis()->GetXmin();
  double shapeEmax = m_Nv->GetYaxis()->GetXmax();
 
  int phbin = m_Nv->GetXaxis()->GetNbins();   
  int ebin = m_Nv->GetYaxis()->GetNbins();

  m_Nv->GetXaxis()->SetLimits(0,m_period);


  //m_Nv->Scale(1.0e+4);  // [ph/(m² s keV)]
  m_flux*= 1.0e+4;

  TH2D *nph = Nph(m_Nv); //ph/m²
  
  int ei2 = nph->GetYaxis()->FindBin(cst::EnNormMin);
  int ei3 = nph->GetYaxis()->FindBin(cst::EnNormMax);

  char temp[200];
  sprintf(temp,"Pulsar shape defined between %.2e  keV and %.2e keV", shapeEmin,shapeEmax);
  WriteToLog(std::string(temp));
  sprintf(temp,"Phase bins: %d ; energy bins: %d",phbin,ebin);
  WriteToLog(std::string(temp));
  sprintf(temp,"Normalization between %.2e keV (%d) and %.2e keV (%d)",cst::EnNormMin,ei2,cst::EnNormMax,ei3);
  WriteToLog(std::string(temp));

  //Normalisation factor according to band between EGRET1 and EGRET2 energies
  //Integration is on a averaged flux over period

  if (NormalizeFlux == 1)
    {
      double norm = m_Nv->Integral(0,phbin,ei2,ei3,"width")/m_period; // ph/m2/s
      m_Nv->Scale(m_flux/norm);
    }



  delete nph;


  //Save output if flags are enabled
  if (::getenv("PULSAR_OUT_NV"))
    {
      SaveNv(m_Nv); // ph/m2/s/keV
    }

  if (SAVETIMEPROFILE)
    {
      SaveTimeProfile(m_Nv); //ph/m2/s/kev
    }

  return m_Nv;
  delete m_Nv;


}


//////////////////////////////////////////////////
TH2D *PulsarSim::Nph(const TH2D *Nv)
{
  TH2D *Nph = (TH2D*) Nv->Clone(); // 1/kev/s
  Nph->SetName("Nph");
  double dei;
  double deltat = Nv->GetXaxis()->GetBinWidth(0);
  
  for (int ei = 0; ei<Ebin; ei++)
    {
      dei   = Nv->GetYaxis()->GetBinWidth(ei+1);
      for(int ti = 0; ti<m_Tbin; ti++)
	{
	  Nph->SetBinContent(ti+1, ei+1, 
			     Nph->GetBinContent(ti+1, ei+1)*dei*deltat); //[1]
	}   
    }
  return Nph;
}

//////////////////////////////////////////////////
/*!
 * \param Nv TH2D ROOT histogram ph/m2/s/keV
 *
 * 
 * This method saves a ROOT file containing the TH2D histogram. The name of the file is created 
 * according to the label of the simulated pulsar.E.g. if the pulsar's name is <i>TEST</i> the output file
 * name is <i>TESTroot.root</i>.If a pre-existent fils exists it will be overwritten. 
*/
void PulsarSim::SaveNv(TH2D *Nv)
{
  
  Nv->SetXTitle("Time [s]");
  Nv->SetYTitle("Energy [keV]");
  Nv->SetZTitle("N_{v} [ph/m^2/s/KeV]");
  Nv->GetXaxis()->SetTitleOffset(1.5);
  Nv->GetYaxis()->SetTitleOffset(1.5);
  Nv->GetZaxis()->SetTitleOffset(1.2);
  Nv->GetXaxis()->CenterTitle();
  Nv->GetYaxis()->CenterTitle();
  Nv->GetZaxis()->CenterTitle();
  
  //Redirect output to a subdirectory

  //Redirect output to a subdirectory
  const char * pulsarOutDir = ::getenv("PULSAROUTFILES");

  // override obssim if running in Gleam environment
  std::string root_name;
  if( pulsarOutDir!=0) 
    root_name = facilities::commonUtilities::joinPath(std::string(pulsarOutDir),  m_name + "root.root");
  else
    root_name =  m_name + "root.root";

  //  std::string root_name = "PsrOutput/" + m_name + "root.root";


  
  TFile mod(root_name.c_str(),"RECREATE");
  Nv->Write();
  mod.Close();
  
}


//////////////////////////////////////////////////
/*!
 * \param Nv TH2D ROOT histogram ph/m2/s/keV
 *
 * 
 * This method saves a Txt file containing the time profile (the TH2D integrated between the max and min energy)
 * The name of the file is created according to the label of the simulated pulsar.
 * E.g. if the pulsar's name is <i>TEST</i> the output file name is <i>TESTTimeProfile.txt</i>. 
*/
void PulsarSim::SaveTimeProfile(TH2D *Nv)
{

  //Redirect output to a subdirectory
  const char * pulsarOutDir = ::getenv("PULSAROUTFILES");

  // override obssim if running in Gleam environment
  std::string nameProfile;
  if( pulsarOutDir!=0) 
    nameProfile = facilities::commonUtilities::joinPath(std::string(pulsarOutDir), m_name + "TimeProfile.txt");
  else
    nameProfile = m_name + "TimeProfile.txt";

  //  std::string nameProfile = "PsrOutput/" + m_name + "TimeProfile.txt";






  ofstream OutTimeProf(nameProfile.c_str());

  int ei2 = Nv->GetYaxis()->FindBin(m_enphmin);
  int ei3 = Nv->GetYaxis()->FindBin(m_enphmax);

  for (int t=0; t < m_Tbin; t++)
    {
      OutTimeProf << t+1 << "\t" <<  Nv->Integral(t+1,t+1,ei2,ei3)*1e4 << "\n";
    }

  OutTimeProf.close();
 
}


//////////////////////////////////////////////////
/*!
 * \param Line line to be written in the output log file
 *
 * This method writes a line into a log file
*/
void PulsarSim::WriteToLog(std::string Line)
{
  std::ofstream PulsarLog(m_LogFileName.c_str(),std::ios::app);
  PulsarLog << "** PulsarSim: " << Line << std::endl;
  PulsarLog.close();
}

