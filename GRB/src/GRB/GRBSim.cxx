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
#define QG false

#define DEBUG 0

using namespace cst;


//////////////////////////////////////////////////
GRBSim::GRBSim(Parameters *params)
  : m_params(params)
{
  m_GRBengine = new GRBengine(params);
  m_fluence   = m_params->GetFluence(); 
  if(DEBUG) 
    std::cout<<" Fluence = "<<m_fluence<<std::endl;
}

void GRBSim::GetUniqueName(const void *ptr, std::string & name)
{
  std::ostringstream my_name;
  my_name << reinterpret_cast<int> (ptr);
  name = my_name.str();
  gDirectory->Delete(name.c_str());

}


TH2D* GRBSim::Fireball()
{
  double *e = new double[Ebin +1];
  for(int i = 0; i<=Ebin; i++)
    {
      e[i] = emin*pow(de,1.0*i); //keV
    }
  
  //////////////////////////////////////////////////  
  
  std::vector<GRBShock*> Shocks = m_GRBengine->CreateShocksVector();
  std::vector<GRBShock*>::iterator pos;
  int nshocks = (int) Shocks.size();
  if(DEBUG)
    {
      std::cout<<"Shock vector size: "<<nshocks<<std::endl;
      for (pos=Shocks.begin();pos!=Shocks.end(); pos++)
	{
	  (*pos)->Print();
	}
    }
  double meanDuration = 0;
  for (pos=Shocks.begin(); pos!=Shocks.end(); pos++)
    {
      (*pos)->SetICComponent(m_params->GetInverseCompton());
      if(DEBUG) (*pos)->Print();
      meanDuration+=(*pos)->GetDuration();
    }
  meanDuration/=nshocks;
  
  //////////////////////////////////////////////////
  double dtqg=0.0;
  double max_tqg=0.0;
  if(QG)
    {
      double dist = GetDistance();
      double Ep = 1e19 * 1e6; //keV;
      dtqg = dist/cst::c/Ep;
      max_tqg = dtqg * 1e7;
    }

  //  int i=1;
  static const double shift =0.0;// Shocks.front()->GetTime() + 2.0*meanDuration;//3.*Shocks.front()->GetDuration();
  m_tfinal = Shocks.back()->GetTime() + meanDuration + shift + max_tqg;

  Tbin     = TMath::Max(100,int(m_tfinal/MinDT));
  if(!m_params->GenerateGBM())
    {
      Tbin     = TMath::Max(100,int(m_tfinal/MinDT));
      Tbin     = TMath::Min(5000,Tbin);  
    }
  if(DEBUG)  
    std::cout<<"Tbin = "<<Tbin<<std::endl;
  gDirectory->Delete("Nv");
  m_Nv = new TH2D("Nv","Nv",Tbin,0.,m_tfinal,Ebin, e);
  double dt = m_Nv->GetBinWidth(1);
  std::string name;
  GetUniqueName(m_Nv,name);
  m_Nv->SetName(name.c_str());
  
  double t = 0.0;
  double energy,tqg;
  
  for(int ti = 0; ti<Tbin; ti++)
    {
      t = ti*dt;
      for(int ei = 0; ei < Ebin; ei++)
	{
	  double nv = m_Nv->GetBinContent(ti+1, ei+1);
	  energy = m_params->rnd->Uniform(e[ei],e[ei+1]);
	  tqg    = t-shift-energy*dtqg;
	  
	  for (int i = 0; i< nshocks; i++)
	    {
	      GRBShock *ashock = Shocks[i];
	      //energy = m_params->rnd->Uniform(e[ei],e[ei+1]);
	      //tqg = t-shift-energy*dtqg
	      nv += ashock->ComputeFlux(tqg,energy);
	    }
	  m_Nv->SetBinContent(ti+1, ei+1, nv);
	  // [ph/(cm² s keV)]
	}
    }
  
  TH2D *nph = Nph(m_Nv); //ph/cm²
  
  int ei1 = nph->GetYaxis()->FindBin(BATSE1);
  int ei2 = nph->GetYaxis()->FindBin(BATSE5);
  double F=0.0;
  double en;
  for (int ei = ei1; ei<=ei2; ei++)
    {
      en   = nph->GetYaxis()->GetBinCenter(ei);
      for(int ti = 1; ti<=Tbin; ti++)
	{
	  F+= nph->GetBinContent(ti, ei) * en;//[keV/cm²]
	}  
    }
  double norm = F*1.0e-3/(erg2meV); //erg/cm²  
  /////////////////////////////////////////////////////////////
  //         IMPORTANT m_Nv has to  be in [ph/(m² s keV)]    //
  //////////////////////////////////////////////////
  m_Nv->Scale(1.0e4 * m_fluence/norm);
  if(DEBUG) 
    m_params->PrintParameters();
  
  delete[] e;
  delete nph;
 
  return m_Nv;
}
//////////////////////////////////////////////////
TH2D *GRBSim::Nph(const TH2D *Nv)
{
  TH2D *Nph = (TH2D*) Nv->Clone(); // [ph/(m² s keV)]  
  std::string name;
  GetUniqueName(Nph,name);
  Nph->SetName(name.c_str());

  double dei;
  double deltat = Nv->GetXaxis()->GetBinWidth(1);
  
  for (int ei = 0; ei<Ebin; ei++)
    {
      dei   = Nv->GetYaxis()->GetBinWidth(ei+1);
      for(int ti = 0; ti<Tbin; ti++)
	{
	  Nph->SetBinContent(ti+1, ei+1, 
			     Nph->GetBinContent(ti+1, ei+1)*dei*deltat); // [ph/(cm²)]  
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
  std::cout<<" Saving "<<root_name<<std::endl;
  TFile mod(root_name,"RECREATE");
  std::string name = m_Nv->GetName();
  m_Nv->SetName("Nv"); // I need a default name.
  m_Nv->Write();
  mod.Close();
  m_Nv->SetName(name.c_str());
}

//////////////////////////////////////////////////

void GRBSim::SaveGBMDefinition(TString GRBname, double ra, double dec, double theta, double phi, double tstart)
{
  TString name = "GRB_";
  name+=GRBname; 
  name+="_DEF.txt";
  std::ofstream os(name,std::ios::out);
  os<<"BURST DEFINITION FILE"<<std::endl;
  os<<"Burst Name"<<std::endl;
  os<<GRBname<<std::endl;
  os<<"RA,DEC (deg):"<<std::endl;
  os<<ra<<" "<<dec<<std::endl;
  os<<"S/C azimuth, elevation (deg):"<<std::endl;
  os<<phi<<" "<<theta<<std::endl;
  os<<"Trigger Time (s):"<<std::endl;
  os<<tstart<<std::endl;
  os.close();
}

double LogFB(double *var, double *par)
{
  
  // GRB function from Band et al.(1993) ApJ.,413:281-292
  double a  = par[0];
  double b  = par[0]+par[1];
  
  double LogE0     = par[2];
  double LogNT     = par[3];
  double LogE      = var[0];
  static const double loge= log10(exp(1.));
  double LogH,LogC; 
  if((a-b)<=0) std::cout<<"WARNING"<<std::endl;
  
  LogH   = log10(a-b) + LogE0;  
  LogC   = (a-b) * (LogH-2.0)-loge*pow(10.0,LogH-LogE0);

  if(LogE <= LogH) 
    return      LogNT + a * (LogE-2.0) - pow(10.0,LogE-LogE0)*loge; 
  return LogC + LogNT + b * (LogE-2.0); // cm^(-2) s^(-1) keV^(-1) 
}



void GRBSim::GetGBMFlux(TString GRBname)
{
  //////////////////////////////////////////////////
  // GBM Spectrum:
  TF1 band("grb_f",LogFB,log10(emin), 4.0, 4); 
  band.SetParNames("a","b","logE0","Log10(Const)");
  
  band.SetParLimits(0,-2.0 , 2.0); // a
  band.SetParLimits(1,-3.0 , -0.001); // b-a; b < a -> b-a < 0 !
  band.SetParLimits(2,log10(emin),4.0);
  
  double a,b,E0,Const;
  TH1D GBM("GBM","GBM",Ebin,log10(emin),log10(emax));
  GBM.SetMinimum(5);
  GBM.SetMaximum(-5);
  
  double t    = 0;
  double dt   = m_Nv->GetXaxis()->GetBinWidth(1);
  double tbin = m_Nv->GetXaxis()->GetNbins();
  
  TString name = "GRB_";
  name+=GRBname; 
  name+="_GBM.txt";
  std::ofstream os(name,std::ios::out);

  os<<"Sample Spectrum File "<<std::endl;
  os<<tbin<<" bins"<<std::endl;
  os<<"Norm   alf   beta  E_p "<<std::endl;
#ifndef WIN32 // THB: Avoid need to link TCanvas on windows 
  if(DEBUG)
    {
      TCanvas *GBMCanvas;
      GBMCanvas = new TCanvas("GBMCanvas","GBMCanvas",500,400);
      std::cout<<"Norm   alf   beta  E_p "<<std::endl;
      GBM.Draw();
    }
#endif
  a =  -1.00;
  b =  -2.25;
  
  
  for(int ti = 0; ti<tbin; ti++)
    {
      t = ti*dt;
      for(int ei = 0; ei < Ebin; ei++)
	{
	  //	  double en = pow(10.0,GBM.GetBinCenter(ei+1));
	  // Notice that Nv is in [ph/(m² s keV)]; nv is in [(ph)/(cm² s keV)]
	  double nv = TMath::Max(1e-10,m_Nv->GetBinContent(ti+1,ei+1)); // [ph/( m² s keV)]
	  nv = log10(nv)-4.0;                                           // [ph/(cm² s keV)]
	  GBM.SetBinContent(ei+1 , nv);                                 // [ph/(cm² s keV)]
	  GBM.SetBinError(ei+1 , nv/100.0);                             // arbitrary small error (1%)
	}
      
      double LogC0  = GBM.GetBinContent(GBM.FindBin(2.0));
      double LogEp0 = 2.0;

      band.SetParameters(a,b-a,LogEp0,LogC0);
      if(LogC0>-5) 
	{
	  if(DEBUG)
	    GBM.Fit("grb_f","r");	  
	  else 
	    GBM.Fit("grb_f","nqr");	  
	} 
      else 
	{
	  band.SetParameters(-2.0,-3.0,LogEp0,-10.0);
	}
      
      a  = band.GetParameter(0);
      b  = a+band.GetParameter(1);
      E0    = pow(10.,band.GetParameter(2));
      Const = pow(10.,band.GetParameter(3));
      double Ep=(2.0+a)*E0; 
      os<<Const<<" "<<a<<" "<<b<<" "<<Ep<<" "<<std::endl;
      
      if(DEBUG)
	{
	  std::cout<<"t= "<<t<<" C= "<<Const<<" a= "<<a<<" b= "<<b<<" E0= "<<E0<<" Ep= "<<Ep<<std::endl;
	  //gPad->SetLogx();
	  //gPad->SetLogy();
	  gPad->Update();
	  //	  TString gbmFlux= "GBMFlux";
	  //	  gbmFlux+=ti;
	  //	  gbmFlux+=".gif";
	  //	  if(ti%10==0) gPad->Print(gbmFlux);
	}
    }   
  os.close();
  //////////////////////////////////////////////////
  //  delete[] e;
}
