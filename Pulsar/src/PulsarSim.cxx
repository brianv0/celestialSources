//////////////////////////////////////////////////
// File PulsarSim.cxx
// contains the code for the implementation of the models
//
//////////////////////////////////////////////////

#include <fstream>
#include <iostream>
#include <ctime>
#include "Pulsar/PulsarConstants.h"
#include "Pulsar/PulsarSim.h"
#include "TFile.h"
#include "TF1.h"

using namespace cst;

//////////////////////////////////////////////////
/*!
 * \param name Name of the pulsar 
 * \param seed Seed of the random generator
 * \param flux Flux of the pulsar in ph/cm2/s 
 * \param enphmin Minimun energy of the extracted photons
 * \param enphmax Maximum energy of the extracted photons
 * \param period Period of the pulsar
 * \param numpeaks Number of peaks
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
 * \param par1 Parameter E0 expressed in GeV
 * \param par2 Parameter En expressed in GeV
 * \param par3 Parameter G 
 * \param par4 Parameter b
 *
 * This method creates a ROOT TH2D histogram according to a phenomenological model. The 2d hist is 
 * obtained by multiplying the lightcurve and the spectrum
 */
TH2D* PulsarSim::PSRPhenom(double par1, double par2, double par3, double par4)
{

  // PSR parameterization for the PolarCap model from Nel & De Jager,Astr&SpSc.230:299-306
  // PSR parameter (Vela)from De Jager et al. 2002
  // The spectrum is generated according to these parameters
  // The lightCurve is random generated with 1 or 2 random peaks whose separation should be greater 
  // than the value of mindist . The amplitudes are random
  // the choice 3 or 4 corresponds to single or double delta fuction.


  // Part 1 - Spectrum

  double En = par1;
  double G1 = par3;
  double E0 = par2;
  double b =  par4;
  double K1 = 138e-8; //This is the constant, will be overwritten by normalization

  std::cout << "\n******** Pulsar Phenomenological Model ********" << std::endl;
  std::cout << "**  Random seed for the model : " << m_seed << std::endl;
  std::cout << "**  Spectrum parameters: " << std::endl;
  std::cout << "**           En = " << En  
	    << " | E0 = " << E0 << std::endl;
  std::cout << "**           G1 = " << G1 
	    << " | b  = "  << b << std::endl;
  std::cout << "**  enphmin " << m_enphmin << " enphmax " << m_enphmax << std::endl; 

  double LowEnBound = TMath::Min(cst::EnNormMin,m_enphmin); 
  if (m_enphmax < cst::EnNormMax)
    m_enphmax = cst::EnNormMax;
  double HighEnBound = TMath::Max(cst::EnNormMax,m_enphmax); 
  
  std::cout << "**           Normalisation between " << cst::EnNormMin << " keV and " << cst::EnNormMax << " keV " << std::endl;
  std::cout << "**           Photon extraction between " << m_enphmin << " keV and " << m_enphmax << " keV " << std::endl;
  std::cout << "**  Spectrum calculated between " << LowEnBound << " keV and " << HighEnBound << " keV " << std::endl; 

  double de = pow(HighEnBound/LowEnBound,1.0/Ebin);
 
  TF1 PulsarSpectralShape("PulsarSpectralShape", 
			  "([0]*((x/[1])^[2])*exp(-1.0*((x/[3])^[4])))", LowEnBound, HighEnBound);

  PulsarSpectralShape.SetParameters(K1,En,G1,E0,b);

  // Part 2- LightCurve

  double fwhm1,fwhm2,peak1,peak2,ampl1,ampl2,mindist;
  mindist = 0.5*m_period; //minimum distance in seconds (minPhase * m_period)
  
  //Set the Random engine.

  TRandom engine;
  //  engine.SetSeed(time(NULL));
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


  //Case of delta Function..

  TH1D deltaFunction("deltaFunction","DeltaFuction",Tbin,0,m_period);

  double dt = m_period/(Tbin-1);

  if ((m_numpeaks == 3 ) || (m_numpeaks == 4 ))
    {
      //      for (int b=0; b < Tbin; b++)
      //{
      //  deltaFunction.SetBinContent(b+1, (ampl1+ampl2)/100);
      //}
      deltaFunction.SetBinContent(int(peak1/dt),ampl1);
      deltaFunction.SetBinContent(int(peak2/dt),ampl2);
    } 


  std::cout << "**\n**  Lightcurve parameters: (midist = " << mindist  << " s.)" << std::endl;
  if (ampl1 !=0)
    {
      std::cout << "**           Peak 1 t = " << peak1 << "(Ph.= " << peak1/m_period << " ) " 
		<< " , FWHM " << fwhm1 << " ampl1 " << ampl1 << std::endl; 
    }
  if (ampl2 !=0)
    {
      std::cout << "**           Peak 2 t = " << peak2 << "(Ph.= " << peak2/m_period << " ) " 
		<< " , FWHM " << fwhm2 << " ampl2 " << ampl2 << std::endl; 
      std::cout << "***********************************************" << std::endl;
    }
  

  TF1 PulsarTimeCurve("PulsarTimeCurve",
		      "([2]*(1/(((x-[0])^2)+(([1]/2)^2))) + [5]*(1/(((x-[3])^2)+(([4]/2)^2))))",
		      0, m_period);
  PulsarTimeCurve.SetParameters(peak1,fwhm1,ampl1,peak2,fwhm2,ampl2);  

  // test in case of stationary curve for flux calibration

  //  TF1 PulsarTimeCurve("PulsarTimeCurve","15.0"); //In order to have a stationary source


  // Part 3 - Combination of Spectrum and lightCurve and filling of TH2D

  double *e = new double[Ebin +1];
  for(int i = 0; i<=Ebin; i++)
    {
      e[i] = LowEnBound*pow(de,1.0*i); //KeV
    }

  gDirectory->Delete("Nv");
  m_Nv = new TH2D("Nv","Nv",Tbin,0.,m_period,Ebin, e);

  //Filling the TH2D Histogram...
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
 
  //Normalisation factor according to band betwee 100MeV and 30 GeV
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

///////////////////////////////////////////////
void PulsarSim::SaveTimeProfile(TH2D *Nv)
{
  // Saving TimeProfile on a TXT Output file.

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


