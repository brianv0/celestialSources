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
  m_fluence                = m_params->GetFluence(); 

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
  while(Shocks[nshocks-i]->GetEfficiency()<1.0e-3)
    { 
      i++;
    }
  m_tfinal = 1.5 * (Shocks[nshocks-i]->GetTime()+Shocks[nshocks-i]->GetDuration());
  double dt = m_tfinal/(Tbin-1);
  gDirectory->Delete("Nv");
  m_Nv = new TH2D("Nv","Nv",Tbin,0.,m_tfinal,Ebin, e);
  
  double t = 0.0;
  
  for (int i = 0; i< nshocks; i++)
    {
      Shocks[i]->SetICComponent(m_params->GetInverseCompton());
      Shocks[i]->Print();
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
	      nv += ashock->ComputeFlux(t,e[ei]);
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
  int ei2 = nph->GetYaxis()->FindBin(BATSE2);
  double norm = nph->Integral(0,Tbin,ei1,ei2,"width")*(1.0e-3)/(erg2meV)/dt; //erg/m²
  
  // IMPORTANT m_Nv has to  be in [ph/(m² s keV)]
  m_Nv->Scale(m_fluence/norm);
  
   
  /*
    for(int i =0; i<(int) Shocks.size();i++)
    delete Shocks[i]; 
  */
  delete e;
  delete nph;
  SaveNv(m_Nv);
  return m_Nv;
}
//////////////////////////////////////////////////
TH2D *GRBSim::Nph(const TH2D *Nv)
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
  Nv->Write();
  mod.Close();
  
  std::ofstream os("grb_generated.txt",std::ios::app);
  os<<m_params->GetGRBNumber()<<" "<<Tmax()<<" "<<m_fluence<<" "<<GRBdir().first<<" "<<GRBdir().second<<std::endl;
  
};

