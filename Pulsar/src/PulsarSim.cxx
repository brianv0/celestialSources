/////////
// Da fare :
// - Controllare che il flusso sia uguale a 0 a t=0 e a t=periodo
// - numero di picchi
//Inserire parametri al modello


#include <fstream>
#include <iostream>
#include <ctime>

#include "Pulsar/PulsarConstants.h"
#include "Pulsar/PulsarSim.h"

#include "TFile.h"
#include "TF1.h"

using namespace cst;

//////////////////////////////////////////////////
PulsarSim::PulsarSim(double fluence,double period, int numpeaks)
{
  m_fluence = fluence; //ph/cm2/s
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
  // The timecurve is generated with 2 Lorenz-shaped curves.
  // The parameters aplitude,fwhm and gamma are random, we require only a minimum distance of 0.4 period


  peak1 = engine->Uniform()*m_period; 
  fwhm1 = engine->Uniform()*m_period;

  int n=0;

  while ((peak1 > 0.5*m_period)
	 || (peak1 < (fwhm1*2)))
    {
      n++;
      peak1 = engine->Uniform()*m_period; 
      fwhm1 = engine->Uniform()*m_period;
      //std::cout << "Extracting peak 1 " << peak1 << " fw1 " << fwhm1 << std::endl;
    }    



  peak2 = engine->Uniform()*m_period; 
  fwhm2 = engine->Uniform()*m_period;


  while ((peak2 > (m_period-2*fwhm2))
	 || (peak2 <(peak1+mindist)))
    {
      peak2 = engine->Uniform()*m_period; 
      fwhm2 = engine->Uniform()*m_period;
      //  std::cout << "Extracting peak 2 " << peak2 << " fw2 " << fwhm2 << std::endl;
      //  std::cout << "Diff " << peak2-peak1 << std::endl; 
    }
  
  ampl1 = engine->Uniform();
  while (ampl1 < 0.3)  ampl1 = engine->Uniform();

  ampl2 = engine->Uniform();
  while (ampl2 < 0.3)  ampl2 = engine->Uniform();
  
  if (m_numpeaks == 1)
    {
      if (engine->Uniform() < 0.5) 
	{
	  ampl1 = 0;
	  //std::cout << "Erasing peak 1 " << std::endl;
	}
      else 
	{
	  ampl2 = 0;
	  // std::cout << "Erasing peak 2 " << std::endl;
	}
    }

  std::cout << "Extracting Random Peaks for Pulsar lightcurve with mindist = " 
	    << mindist  << " msec. " << std::endl;
  std::cout << "----> LightCurve " << std::endl;
  std::cout << "      Peak 1 t = " << peak1 << "(Phase= " << peak1/m_period << " ) " 
	    << " FWHM " << fwhm1 << " ampl1 " << ampl1 << std::endl;  
  std::cout << "      Peak 2 t = " << peak2 << "(Phase= " << peak2/m_period << " ) " 
	    <<" FWHM " << fwhm2 << " ampl2 " << ampl2 << std::endl; 

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


  //  Part 2 - Generation of the Spectrum.
  //
  
  
  double K1 = 138e-8; //This is the constant, irrelevant because there is the m_fluence normalizatio
 
  TF1 *PulsarSpectralShape = new TF1("PulsarSpectralShape", 
				     "([0]*((x/[1])^[2])*exp(-1.0*((x/[3])^[4])))",cst::emin,cst::emax);

  PulsarSpectralShape->SetParameters(K1,En,G1,E0,b);



  //Fill the TH2D object

  double t = 0.0;
  for(int ti = 0; ti<Tbin; ti++)
    {
      t = ti*dt;
      double nt = PulsarTimeCurve->Eval(t);
      for(int ei = 0; ei < Ebin; ei++)
	{
	  double nv = PulsarSpectralShape->Eval(e[ei]);
	  m_Nv->SetBinContent(ti+1, ei+1, nt*nv);// [ph/(cm² s KeV)]
	}
    }


  // Conervsion 1/cm² -> 1/m² IMPORTANT m_Nv has to be in [ph/m2s KeV)
  // In the XML file is espressed in Kev/cm2, but here we have to use erg/m2

  
  m_Nv->Scale(1.0e+4);        // [ph/(m² s keV)]
  m_fluence *= 1.0e+4;        //  ph/m²/s
  
  // nph = nv * dE * dt


  TH2D *nph = Nph(m_Nv); //ph/m²
  
  //  int ei1 = nph->GetYaxis()->FindBin(BATSE1);
  //  int ei2 = nph->GetYaxis()->FindBin(BATSE2);
  //double norm = nph->Integral(0,Tbin,ei1,ei2,"width")*(1.0e-3)/(dt*erg2meV); //erg/m²

  //Normalisation with m_fluence
  // m_fluence is now converted in erg/m2, so we must convert the integral into erg/m2
  // IMPORTANT m_Nv has to  be in [ph/(m² s keV)]

  /* now we want to normalize in ph/cm2/s and do the underliying lines are not required
    double norm = nph->Integral(0,Tbin,0,Ebin,"width")*(1.0e-3)/(dt*erg2meV); //erg/m²
  */

  int ei2 = nph->GetYaxis()->FindBin(EGRET2);
  int ei3 = nph->GetYaxis()->FindBin(EGRET3);
 
  double norm = nph->Integral(0,Tbin,ei2,ei3)/m_period;//*(1.0e-3)/(dt*erg2meV); //ph/(m² s)

  m_Nv->Scale(m_fluence/norm);

  delete e;
  delete nph;
  SaveNv(m_Nv);
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

