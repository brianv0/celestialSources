#include <fstream>
#include <iostream>

#include "GRBConstants.h"
#include "GRBShell.h"
#include "GRBShock.h"
#include "GRBengine.h"
#include "GRBSim.h"

#include "TFile.h"

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
  int nshocks = (int) Shocks.size();
  int i=1;
  double shift = 3.*Shocks.front()->GetDuration();
  m_tfinal = 1.5*Shocks.back()->GetTime()+shift;
  std::cout<<shift<<" "<<m_tfinal<<" "<<Shocks.back()->GetTime()<<std::endl;

  double dt = m_tfinal/(Tbin-1);
  
  m_Nv = new TH2D("Nv","Nv",Tbin,0.,m_tfinal,Ebin, e);
  std::string name;
  GetUniqueName(m_Nv,name);
  gDirectory->Delete(name.c_str());
  m_Nv->SetName(name.c_str());
  
  double t = 0.0;
  
  for (int i = 0; i< nshocks; i++)
    {
      Shocks[i]->SetICComponent(m_params->GetInverseCompton());
      //      Shocks[i]->Print();
    }

  
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
void GRBSim::SaveNv(TH2D *Nv)
{
  
  Nv->SetXTitle("Time [s]");
  Nv->SetYTitle("Energy [keV]");
  Nv->SetZTitle("N_{v} [ph/m^2/s/keV]");
  Nv->GetXaxis()->SetTitleOffset(1.5);
  Nv->GetYaxis()->SetTitleOffset(1.5);
  Nv->GetZaxis()->SetTitleOffset(1.2);
  Nv->GetXaxis()->CenterTitle();
  Nv->GetYaxis()->CenterTitle();
  Nv->GetZaxis()->CenterTitle();
  
  char root_name[100];
  sprintf(root_name,"grb_%d.root",(int)m_params->GetGRBNumber());
  TFile mod(root_name,"RECREATE");
  std::string name = Nv->GetName();
  Nv->SetName("Nv"); // I need a default name.
  Nv->Write();
  Nv->SetName(name.c_str());
  mod.Close();
};

