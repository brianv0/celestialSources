#include <fstream>
#include <iostream>

#include "GRBConstants.h"
#include "GRBShell.h"
#include "GRBShock.h"
#include "GRBengine.h"
#include "GRBSim.h"

#include "TFile.h"

#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"

#define DEBUG 0 

using namespace cst;


//////////////////////////////////////////////////
GRBSim::GRBSim(Parameters *params)
  : m_params(params)
{
  m_GRBengine = new GRBengine(params);
  m_fluence   = m_params->GetFluence(); 
}

void GRBSim::GetUniqueName(const void *ptr, std::string & name)
{
  std::ostringstream my_name;
  my_name << reinterpret_cast<int> (ptr);
  name = my_name.str();
}


TH2D* GRBSim::Fireball()
{
  double *e = new double[Ebin +1];
  for(int i = 0; i<=Ebin; i++)
    {
      e[i] = emin*pow(de,1.0*i); //keV
    }
  
  //////////////////////////////////////////////////  
  
  std::vector<GRBShock*> Shocks = 
    m_GRBengine->CreateShocksVector();

  double meanDuration = 0;
  int nshocks = (int) Shocks.size();

  std::vector<GRBShock*>::iterator pos;
  for (pos=Shocks.begin();pos!=Shocks.end();++pos)
    {
      (*pos)->SetICComponent(m_params->GetInverseCompton());
      if(DEBUG) (*pos)->Print();
      meanDuration+=(*pos)->GetDuration();
    }
  meanDuration/=nshocks;

  //  int i=1;
  double shift =0.0;// Shocks.front()->GetTime() + 2.0*meanDuration;//3.*Shocks.front()->GetDuration();
  m_tfinal = Shocks.back()->GetTime() + meanDuration + shift;
  if(DEBUG) std::cout<<shift<<" "<<m_tfinal<<" "<<meanDuration<<std::endl;
  double dt = m_tfinal/(Tbin-1);
  
  m_Nv = new TH2D("Nv","Nv",Tbin,0.,m_tfinal,Ebin, e);
  std::string name;
  GetUniqueName(m_Nv,name);
  gDirectory->Delete(name.c_str());
  m_Nv->SetName(name.c_str());
  
  double t = 0.0;
  
  for(int ti = 0; ti<Tbin; ti++)
    {
      t = ti*dt;
      for(int ei = 0; ei < Ebin; ei++)
	{
	  double nv = m_Nv->GetBinContent(ti+1, ei+1);
	  for (int i = 0; i< nshocks; i++)
	    {
	      GRBShock *ashock = Shocks[i];
	      nv += ashock->ComputeFlux(t-shift,e[ei]);
	    }
	  m_Nv->SetBinContent(ti+1, ei+1, nv);
	  // [ph/(cm² s keV)]
	}
    }
  // Conersion 1/cm² -> 1/m²
  m_Nv->Scale(1.0e+4); // [ph/(m² s keV)]
  // double fluence = 1.0e-6; //erg/cm²
  m_fluence *= 1.0e+4;        //erg/m²
  // nph = nv * dE * dt
  TH2D *nph = Nph(m_Nv); //ph/m²
  
  int ei1 = nph->GetYaxis()->FindBin(BATSE1);
  int ei2 = nph->GetYaxis()->FindBin(BATSE5);
  double norm = nph->Integral(0,Tbin,ei1,ei2,"width")*(1.0e-3)/(erg2meV)/dt; //erg/m²
  
  if (norm<1e-20) 
    {
      m_params->SetInitialSeparation(m_params->GetInitialSeparation()/5.0);
      m_fluence =  m_params->GetBATSEFluence();
      m_params->SetFluence(m_fluence);
      delete[] e;
      return Fireball();
    }
  
  // IMPORTANT m_Nv has to  be in [ph/(m² s keV)]
  
  m_Nv->Scale(m_fluence/norm);
  
   
  Shocks.erase(Shocks.begin(), Shocks.end());
  
  if(DEBUG) 
    m_params->PrintParameters();

  delete[] e;
  delete nph;
  //SaveNv(m_Nv);
  return m_Nv;
}
//////////////////////////////////////////////////
TH2D *GRBSim::Nph(const TH2D *Nv)
{
  TH2D *Nph = (TH2D*) Nv->Clone(); // [ph/(m² s keV)]  
  std::string name;
  GetUniqueName(Nph,name);
  gDirectory->Delete(name.c_str());
  Nph->SetName(name.c_str());

  double dei;
  double deltat = Nv->GetXaxis()->GetBinWidth(1);
  
  for (int ei = 0; ei<Ebin; ei++)
    {
      dei   = Nv->GetYaxis()->GetBinWidth(ei+1);
      for(int ti = 0; ti<Tbin; ti++)
	{
	  Nph->SetBinContent(ti+1, ei+1, 
			     Nph->GetBinContent(ti+1, ei+1)*dei*deltat); //[ph]
	}   
    }
  return Nph;
}

//////////////////////////////////////////////////
void GRBSim::SaveNv()
{
  
  m_Nv->SetXTitle("Time [s]");
  m_Nv->SetYTitle("Energy [keV]");
  m_Nv->SetZTitle("N_{v} [ph/m^2/s/keV]");
  m_Nv->GetXaxis()->SetTitleOffset(1.5);
  m_Nv->GetYaxis()->SetTitleOffset(1.5);
  m_Nv->GetZaxis()->SetTitleOffset(1.2);
  m_Nv->GetXaxis()->CenterTitle();
  m_Nv->GetYaxis()->CenterTitle();
  m_Nv->GetZaxis()->CenterTitle();
  
  char root_name[100];
  sprintf(root_name,"grb_%d.root",(int)m_params->GetGRBNumber());
  
  TFile mod(root_name,"RECREATE");
  std::string name = m_Nv->GetName();
  m_Nv->SetName("Nv"); // I need a default name.
  m_Nv->Write();
  m_Nv->SetName(name.c_str());
  mod.Close();
};

//////////////////////////////////////////////////
double BandF(double *var, double *par)
{
  double a  = par[0];
  double b  = par[1];
  double E0 = pow(10.,par[2]);
  double NT = pow(10.,par[3]);
  
  double E   = var[0];
  double C   = pow((a-b)*E0,a-b)*exp(b-a);
  double H   = (a-b) * E0;
  if(E <= H) 
    return NT *E* E* pow(E,a) * exp(-E/E0);
  return C* NT *E* E* pow(E,b); // ph cm^(-2) s^(-1) keV
}

void GRBSim::GetGBMFlux()
{
  double *e = new double[Ebin +1];
  for(int i = 0; i<=Ebin; i++)
    {
      e[i] = emin*pow(de,1.0*i); //keV
    }
  //////////////////////////////////////////////////
  // GBM Spectrum:
  TF1 band("grb_f",BandF,emin,1.0e+4,4); 
  band.SetParNames("a","b","Log10(E0)","Log10(Const)");
  band.SetParameters(-1.4,-2.4,2.5,-3.0);
  band.SetParLimits(0,-9.0,0.0);
  band.SetParLimits(1,-10.0,0.0);
  band.SetParLimits(2,log10(emin),4.0);
  band.SetParLimits(3,-5.0,5.0);
  double a,b,E0,Const;
  TH1D GBM("GBM","GBM",Ebin, e);
  GBM.SetMinimum(1e-5);
  GBM.SetMaximum(1e6);
  double t=0;
  double dt = m_Nv->GetXaxis()->GetBinWidth(1);
  for(int ti = 0; ti<Tbin; ti++)
    {
      t = ti*dt;
      for(int ei = 0; ei < Ebin; ei++)
	{
	  double nv = m_Nv->GetBinContent(ti+1, ei+1); // [ph/(m² s keV)]
	  GBM.SetBinContent(ei+1, e[ei]*e[ei]* nv *1.0e-4); // [(keV*keV)/(cm² s keV)]
	}
      GBM.Fit("grb_f","rq");
      
      a=band.GetParameter(0);
      b=band.GetParameter(1);
      E0=pow(10.,band.GetParameter(2));
      Const=pow(10.,band.GetParameter(3));
      band.SetParameters(a,b,band.GetParameter(2),band.GetParameter(3));
      //Ep=(a+2)*E0;
      std::cout<<t<<" "<<a<<" "<<b<<" "<<E0<<" "<<Const<<std::endl;
      gPad->SetLogx();
      gPad->SetLogy();
      gPad->Update();
      TString gbmFlux= "GBMFlux";
      gbmFlux+=ti;
      gbmFlux+=".gif";
      if(ti%10==0) gPad->Print(gbmFlux);
    }
  //////////////////////////////////////////////////
  delete[] e;
}
