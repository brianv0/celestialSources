//////////////////////////////////////////////////
// File PulsarSim.cxx
// contains the code for the implementation of the models
//////////////////////////////////////////////////

#include <fstream>
#include <iostream>
#include <iomanip>
#include <ctime>
#include "Pulsar/PulsarConstants.h"
#include "Pulsar/PulsarSim.h"
#include "TFile.h"
#include "TF1.h"

using namespace cst;

//////////////////////////////////////////////////
/*!
 * \param name Name of the pulsar; 
 * \param seed Seed of the random generator;
 * \param flux Flux of the pulsar in ph/cm2/s; 
 * \param enphmin Minimun energy of the extracted photons in keV;
 * \param enphmax Maximum energy of the extracted photons in keV;
 * \param period Period of the pulsar;
 * \param numpeaks Number of peaks : 1 - Only one peak;
 *                                   2 - Two peaks;
 *                                   3 - One delta shaped peak (for test purposes);
 *                                   4 - Two delta shaped peaks (for test purposes);
 * 
 */
PulsarSim::PulsarSim(std::string name, int seed, double flux, double enphmin, double enphmax, double period, int numpeaks)
{

  m_flux = flux; //ph/cm2/s
  m_period  = period;
  m_numpeaks = numpeaks;
  m_enphmin = enphmin;
  m_enphmax = enphmax;
  m_name = name;
  m_seed = seed;
}

//////////////////////////////////////////////////
/*!
 * \param par1 Parameter E0 expressed in GeV;
 * \param par2 Parameter En expressed in GeV;
 * \param par3 Parameter g ;
 * \param par4 Parameter b ;
 *
 * This method creates a ROOT TH2D histogram according to a phenomenological model. The 2d hist is 
 * obtained by multiplying the lightcurve and the spectrum. The lightcurve is obtained by random generating 
 * a profile of 1 or 2 peaks separated by a minimum distance, that currently is set to one half of the period. 
 * The spectrum is generated from an analytical form :
 *
 * \image html NJSnForm.gif 
 * <br>
 * where the parameters E0,En,g and g are the parameters of the method. This form describes a power law spectrum
 * with an exponential cutoff. This formula is choosen according to Nel, De Jager (1995,see Ref.below ) and 
 * De Jager (2002, see Ref.below). According to Cheng (1994, see Ref. below) an Outer Gap scenario can be obtained 
 * by setting b=1.
 * The constant K indicates the normalisation. We choose to set the normalisation in accord to the 3rd EGRET Catalog
 * where the fluxes are reported above 100MeV in ph/s/cm2. In our simulator we adopt the same convention.
 * The ROOT histogram is then saved to a ROOT file with the same name of the pulsar, and also a Txt Time profile id
 * saved.
 * For more informations and fot a brief tutorial please see:
 * <br>
 * <a href="#dejager02">http://www.pi.infn.it/~razzano/Pulsar/PulsarSpectrumTutorial/PsrSpectrumTut.html </a>
 * <br><br>
 * <i><b>References:</b></i>
 * <ul>
 *  <li><a name="chengpaper94"></a>Cheng, K.S. and Ding, W.K.Y.:1994, <i>Astrophys.J.</i>431,724;</li>
 *  <li><a name="neilpaper95"></a>Neil, H.I. and De Jager, O.C.:1995, <i>Astrophysics and Space Science</i> 230:209-306;</li>
 *  <li><a name="dejager02"></a>De Jager, O.C.:2002, <i>African Skies</i>,No 7;</li>
 *  <li><a name="egret3cat"></a>Hartman, R.C. et al.: 1999, <i>The Astrophysical Journal Supplement Series</i>,123:79-202;</li>
 * </ul>
 */
TH2D* PulsarSim::PSRPhenom(double par1, double par2, double par3, double par4)
{

  // Part 1 - Spectrum

  double En = par1;
  double G1 = par3;
  double E0 = par2;
  double b =  par4;
  double K1 = 138e-8; //This is the constant, will be overwritten by normalization

  //Write out informations about the model parameters
  std::cout << "\n******** Pulsar Phenomenological Model ********" << std::endl;
  std::cout << "**  Random seed for the model : " << m_seed << std::endl;
  std::cout << "**  Spectrum parameters: " << std::endl;
  std::cout << "**           En = " << En  
	    << " | E0 = " << E0 << std::endl;
  std::cout << "**           G1 = " << G1 
	    << " | b  = "  << b << std::endl;
  std::cout << "**  enphmin " << m_enphmin << " enphmax " << m_enphmax << std::endl; 

  //Establish the lower and upper limits for the energy in ROOT histogram.Write out these infos
  double LowEnBound = TMath::Min(cst::EnNormMin,m_enphmin); 
  if (m_enphmax < cst::EnNormMax)
    m_enphmax = cst::EnNormMax;
  double HighEnBound = TMath::Max(cst::EnNormMax,m_enphmax); 
  
  std::cout << "**           Normalisation between " << cst::EnNormMin << " keV and " 
	    << cst::EnNormMax << " keV " << std::endl;
  std::cout << "**           Photon extraction between " << m_enphmin << " keV and " 
	    << m_enphmax << " keV " << std::endl;
  std::cout << "**  Spectrum calculated between " << LowEnBound << " keV and " 
	    << HighEnBound << " keV " << std::endl; 

  //Create the spetrum profile
  double de = pow(HighEnBound/LowEnBound,1.0/Ebin);
 
  TF1 PulsarSpectralShape("PulsarSpectralShape", 
			  "([0]*((x/[1])^[2])*exp(-1.0*((x/[3])^[4])))", LowEnBound, HighEnBound);
  PulsarSpectralShape.SetParameters(K1,En,G1,E0,b);



  // Part 2- LightCurve

  double fwhm1,fwhm2,peak1,peak2,ampl1,ampl2,mindist;
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
  if ((m_numpeaks == 1) || (m_numpeaks == 3))
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

  //Case of delta Function (numPeaks =3 or =4)
  TH1D deltaFunction("deltaFunction","DeltaFuction",Tbin,0,m_period);

  double dt = m_period/(Tbin-1);

  if ((m_numpeaks == 3 ) || (m_numpeaks == 4 ))
    {
      deltaFunction.SetBinContent(int(peak1/dt),ampl1);
      deltaFunction.SetBinContent(int(peak2/dt),ampl2);
    } 

  //writes out informations about hte lightcurves.
  std::cout << std::setprecision(3) << "**\n**  Lightcurve parameters: (midist = " 
	    << mindist  << " s.)" << std::endl;
  if (ampl1 !=0)
    {
      std::cout << std::setprecision(3) << "**           Peak 1 t = " << peak1 << "(Ph.= " << peak1/m_period << " ) " 
		<< " , FWHM " << fwhm1 << " ampl1 " << ampl1 << std::endl; 
    }
  if (ampl2 !=0)
    {
      std::cout << std::setprecision(3) << "**           Peak 2 t = " << peak2 << "(Ph.= " << peak2/m_period << " ) " 
		<< " , FWHM " << fwhm2 << " ampl2 " << ampl2 << std::endl; 
      std::cout << "***********************************************" << std::endl;
    }
  

  TF1 PulsarTimeCurve("PulsarTimeCurve",
		      "([2]*(1/(((x-[0])^2)+(([1]/2)^2))) + [5]*(1/(((x-[3])^2)+(([4]/2)^2))))",
		      0, m_period);
  PulsarTimeCurve.SetParameters(peak1,fwhm1,ampl1,peak2,fwhm2,ampl2);  
  
  // test curve in case of stationary curve for flux calibration
  // TF1 PulsarTimeCurve("PulsarTimeCurve","15.0"); //In order to have a stationary source
    
  // Part 3 - Combination of Spectrum and lightCurve and filling of TH2D

  double *e = new double[Ebin +1];
  for(int i = 0; i<=Ebin; i++)
    {
      e[i] = LowEnBound*pow(de,1.0*i); //KeV
    }

  gDirectory->Delete("Nv");
  m_Nv = new TH2D("Nv","Nv",Tbin,0.,m_period,Ebin, e);

  //Filling the TH2D Histogram
  double t = 0.0;
  for(int ti = 0; ti<Tbin; ti++)
    {
      t = ti*dt;
      double nt = PulsarTimeCurve.Eval(t);
      if ((m_numpeaks == 3 ) || (m_numpeaks == 4 ))
	{
	  nt = deltaFunction.GetBinContent(ti);
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
  double norm = m_Nv->Integral(0,Tbin,ei2,ei3,"width")/m_period; // ph/m2/s

  m_Nv->Scale(m_flux/norm);

  delete nph;
  delete[] e;
  
  SaveNv(m_Nv); // ph/m2/s/keV
  SaveTimeProfile(m_Nv); //ph/m2/s/kev

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
      for(int ti = 0; ti<Tbin; ti++)
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
 * according to the label of the simulated pulsar.E.g. if the pulsar's name is <i>Test</i> the output file
 * name is <i>Testroot.root</i>.If a preexistent fils exists it will be overwritten. 
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
  
  char root_name[100];

  for (int i=0; i< m_name.length()+1; i++)
    {
      root_name[i] = m_name[i];
    }

  sprintf(root_name,"%sroot.root",root_name);
  
  TFile mod(root_name,"RECREATE");
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
 * E.g. if the pulsar's name is <i>Test</i> the output file
 * name is <i>TestTimeProfile.txt</i>. 
*/
void PulsarSim::SaveTimeProfile(TH2D *Nv)
{

  char temp[30];

  for (int i=0; i< m_name.length()+1; i++)
    {
      temp[i] = m_name[i];
    }

  sprintf(temp,"%sTimeProfile.txt",temp);
  ofstream OutTimeProf(temp);

  int ei2 = Nv->GetYaxis()->FindBin(m_enphmin);
  int ei3 = Nv->GetYaxis()->FindBin(m_enphmax);

  for (int t=0; t < Tbin; t++)
    {
      OutTimeProf << t+1 << "\t" <<  Nv->Integral(t+1,t+1,ei2,ei3)*1e4 << "\n";
    }

  OutTimeProf.close();

};


