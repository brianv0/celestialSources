/////////
// Da fare :
// - Controllare che il flusso sia uguale a 0 a t=0 e a t=periodo

//edds:eeciao
#include <fstream>
#include <iostream>
#include <ctime>

#include "Pulsar/PulsarConstants.h"
#include "Pulsar/PulsarSim.h"

#include "TFile.h"
#include "TF1.h"

using namespace cst;

//////////////////////////////////////////////////
PulsarSim::PulsarSim(double flux,double period, int numpeaks)
{
  m_flux = flux; //ph/cm2/s
  m_period  = period;
  m_numpeaks = numpeaks;
}

TH2D* PulsarSim::PSRPhenom(double par1, double par2, double par3, double par4)
{

  // PSR parameterization for the PolarCap model from Nel & De Jager,Astr&SpSc.230:299-306
  // PSR parameter (Vela)from De Jager et al. 2002

  double En = par1;
  double G1 = par3;
  double E0 = par2;
  double b =  par4;
  std::cout << "----> Spectral profile parameters " << std::endl;
  std::cout << "      En = " << En  
	    << " | E0 = " << E0 << std::endl;
  std::cout << "      G1 = " << G1 
	    << " | b  = "  << b << std::endl;

  double fwhm1,fwhm2,peak1,peak2,ampl1,ampl2,mindist;
  mindist = 0.5*m_period; //minimum distance in seconds (minPhase * m_period)
  
  //Set the Random engine.

  TRandom *engine = new TRandom();
  //engine->SetSeed(time(NULL));
  engine->SetSeed(0);
  
  // Part 1 - Generate a single or double peak-shaped timecurve with Lorenz function.
  // The timecurve is generated with 1 or 2 Lorenz-shaped curves.
  // The parameters amplitude,fwhm and gamma are random, we require only a minimum distance of 0.4 period
  // If the user select 3 or 4 is possible to have a time profile with 1 or 2 delta instead of the
  // usual Lorenzian profile.


  peak1 = engine->Uniform()*m_period; 
  fwhm1 = engine->Uniform()*m_period;

  int n=0;

  while ((peak1 > 0.5*m_period)
	 || (peak1 < (fwhm1*2)))
    {
      n++;
      peak1 = engine->Uniform()*m_period; 
      fwhm1 = engine->Uniform()*m_period;
    }    



  peak2 = engine->Uniform()*m_period; 
  fwhm2 = engine->Uniform()*m_period;


  while ((peak2 > (m_period-2*fwhm2))
	 || (peak2 <(peak1+mindist)))
    {
      peak2 = engine->Uniform()*m_period; 
      fwhm2 = engine->Uniform()*m_period;
    }
  
  ampl1 = engine->Uniform();
  while (ampl1 < 0.3)  ampl1 = engine->Uniform();

  ampl2 = engine->Uniform();
  while (ampl2 < 0.3)  ampl2 = engine->Uniform();
  
  if ((m_numpeaks == 1) || (m_numpeaks == 3))
    {
      if (engine->Uniform() < 0.5) 
	{
	  ampl1 = 0;
	}
      else 
	{
	  ampl2 = 0;
	}
    }

  std::cout << "\n*****\nExtracting Random LightCurve with mindist = " 
	    << mindist  << " msec. " << std::endl;
  std::cout << "----> LightCurve " << std::endl;
  std::cout << "      Peak 1 t = " << peak1 << "(Phase= " << peak1/m_period << " ) " 
	    << " FWHM " << fwhm1 << " ampl1 " << ampl1 << std::endl;  
  std::cout << "      Peak 2 t = " << peak2 << "(Phase= " << peak2/m_period << " ) " 
	    <<" FWHM " << fwhm2 << " ampl2 " << ampl2 << "\n*****" << std::endl; 

  TF1 *PulsarTimeCurve = new TF1("PulsarTimeCurve",
                                 "([2]*(1/(((x-[0])^2)+(([1]/2)^2))) + [5]*(1/(((x-[3])^2)+(([4]/2)^2))))",
                                 0, m_period);
  PulsarTimeCurve->SetParameters(peak1,fwhm1,ampl1,peak2,fwhm2,ampl2);  


  double *e = new double[Ebin +1];
  for(int i = 0; i<=Ebin; i++)
    {
      e[i] = emin*pow(de,1.0*i); //KeV
    }

  double dt = m_period/(Tbin-1);
  gDirectory->Delete("Nv");
  m_Nv = new TH2D("Nv","Nv",Tbin,0.,m_period,Ebin, e);


  //  Part 2 - Generation of the Spectrum.according to the Phenomenological model
  
  
  double K1 = 138e-8; //This is the constant, irrelevant because there is the m_fluxnormalizatio
 
  TF1 *PulsarSpectralShape = new TF1("PulsarSpectralShape", 
				     "([0]*((x/[1])^[2])*exp(-1.0*((x/[3])^[4])))",cst::emin,cst::emax);

  PulsarSpectralShape->SetParameters(K1,En,G1,E0,b);


  // This deltaFunction is used instead of the usual curve if the user select 3 or 4 
  // as parameters for number of peaks.

  TH1D *deltaFunction = new TH1D("deltaFunction","DeltaFuction",Tbin,0,m_period);

  if ((m_numpeaks == 3 ) || (m_numpeaks == 4 ))
    {
      
      deltaFunction->SetBinContent(int(peak1/dt),ampl1);
      deltaFunction->SetBinContent(int(peak2/dt),ampl2);
      for (int i=0; i < Tbin; i++)
	deltaFunction->SetBinContent(i,deltaFunction->GetBinContent(i)+((ampl1+ampl2)/100));
    } 


  //Filling the TH2D Histogram...
  double t = 0.0;
  for(int ti = 0; ti<Tbin; ti++)
    {
      t = ti*dt;
      double nt = PulsarTimeCurve->Eval(t);

      if ((m_numpeaks == 3 ) || (m_numpeaks == 4 ))
	{
	  nt = deltaFunction->GetBinContent(ti);
	}

      for(int ei = 0; ei < Ebin; ei++)
	{
	  double nv = PulsarSpectralShape->Eval(e[ei]);
	  m_Nv->SetBinContent(ti+1, ei+1, nt*nv);// [ph/(cm² s KeV)]
	}
    }


  // Conversion 1/cm² -> 1/m² IMPORTANT m_Nv has to be in [ph/m2s KeV)
  // In the XML file the flux is espressed in ph/cm2/s according to EGRET Catalogs

  m_Nv->Scale(1.0e+4);  // [ph/(m² s keV)]
  m_flux*= 1.0e+4;      //  ph/m²/s
  
  // nph = nv * dE * dt

  TH2D *nph = Nph(m_Nv); //ph/m²
  
  int ei2 = nph->GetYaxis()->FindBin(EGRET2);
  int ei3 = nph->GetYaxis()->FindBin(EGRET3);
 
  //Normalisation Integration above 100MeV up to limit of EGRET band 
  double norm = nph->Integral(0,Tbin,ei2,ei3)/m_period; //ph/cm2/s

  m_Nv->Scale(m_flux/norm);

  delete e;
  delete nph;
  delete deltaFunction;
  SaveNv(m_Nv);

  //Saving TimeProfile on a TXT Output file.

  ofstream OutTimeProf("PSRTimeProfile.txt");
  OutTimeProf <<  m_flux << "\t" << m_period << "\n-----------------------\n";  
  TH1D *NvTimeProf = new TH1D("NvTimeProf","NvTimeProf",Tbin,0.0,m_period);

  for (int t=0; t < Tbin; t++)
    {
      NvTimeProf->SetBinContent(t+1,m_Nv->Integral(t+1,t+1,ei2,ei3));
      OutTimeProf << t << "\t" <<  NvTimeProf->GetBinContent(t+1)*1e4 << "\n";
    }


  OutTimeProf.close();


  return m_Nv;
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
  sprintf(root_name,"pulsar.root");
  
  TFile mod(root_name,"RECREATE");
  Nv->Write();
  mod.Close();
  
};


