#include "GRBConstants.h"
#include "GRBShell.h"
#include "GRBShock.h"
#include "GRBengine.h"
#include "GRBSim.h"

using namespace cst;

//////////////////////////////////////////////////
GRBSim::GRBSim(Parameters *params)
  : m_params(params)
{
  m_GRB = new GRBengine(params);
}

TH2D* GRBSim::Fireball()
{
  int    Nshell            = m_params->m_nshell;
  double BATSE_fluence     = m_params->m_fluence;
  double Etot              = m_params->m_etot;
  double InitialSeparation = m_params->m_initialSeparation;
  double InitialThickness  = m_params->m_initialThickness;
  
  double *e = new double[Ebin +1];
  for(int i = 0; i<=Ebin; i++)
    {
      e[i] = emin*pow(de,1.0*i); //keV
    }
  
  //////////////////////////////////////////////////  
  
  std::vector<GRBShock*> Shocks = 
    m_GRB->CreateShocksVector(Nshell,InitialSeparation,InitialThickness,Etot);
  int nshocks = (int) Shocks.size();
  //  if(nshocks==0) return;
  
  m_tfinal = 1.5 * Shocks[nshocks-1]->GetTime();
  double dt = m_tfinal/(Tbin-1);
  gDirectory->Delete("Nv");
  m_Nv = new TH2D("Nv","Nv",Tbin,0.,m_tfinal,Ebin, e);
  
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
	      nv += ashock->ComputeFlux(t,e[ei]);
	    }
	  m_Nv->SetBinContent(ti+1, ei+1, nv);
	  // [ph/(cm² s keV)]
	}
    }
  // Conersion 1/cm² -> 1/m²
  m_Nv->Scale(1.0e+4); // [ph/(m² s keV)]
  // double fluence = 1.0e-6; //erg/cm²
  BATSE_fluence *= 1.0e+4;        //erg/m²
  // nph = nv * dE * dt
  TH2D *nph = Nph(m_Nv); //ph/m²
  
  int ei1 = nph->GetYaxis()->FindBin(BATSE1);
  int ei2 = nph->GetYaxis()->FindBin(BATSE2);
  double norm = nph->Integral(0,Tbin,ei1,ei2,"width")*(1.0e-3)/(dt*erg2meV); //erg/m²
  
  // IMPORTANT m_Nv has to  be in [ph/(m² s keV)]
  m_Nv->Scale(BATSE_fluence/norm);  
  /*
    nph = Nph(m_Nv); //ph/m²
    cout<<nph->Integral(0,Tbin,ei1,ei2,"width")*(1.0e-7)/(dt*erg2meV)<<endl; //erg/cm²
  */    
  delete nph;
  SaveNv();
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
  
  TFile *mod = new TFile("grb.root","RECREATE");
  m_Nv->Write();
  mod->Close();
};

