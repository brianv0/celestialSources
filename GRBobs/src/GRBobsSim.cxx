#include <fstream>
#include <iostream>

#include "GRBobs/GRBobsConstants.h"
#include "GRBobs/GRBobsengine.h"
#include "GRBobs/GRBobsSim.h"
#include "GRBobs/GRBobsPulse.h"

#include "TFile.h"
#include "TCanvas.h"
#define DEBUG 0

using namespace ObsCst;
//////////////////////////////////////////////////
GRBobsSim::GRBobsSim(GRBobsParameters *params)
  : m_params(params)
{
  m_GRBengine = new GRBobsengine(params);
}

void GRBobsSim::GetUniqueName(const void *ptr, std::string & name)
{
  std::ostringstream my_name;
  my_name << reinterpret_cast<int> (ptr);
  name = my_name.str();
  gDirectory->Delete(name.c_str());
}

TH2D* GRBobsSim::MakeGRB()
{
  double *e = new double[Ebin +1];
  for(int i = 0; i<=Ebin; i++)
    {
     e[i] = emin*pow(de,1.0*i); //keV
    }
  
  //////////////////////////////////////////////////  
  double s_TimeBinWidth=TimeBinWidth;
  std::vector<GRBobsPulse*> Pulses = m_GRBengine->CreatePulsesVector();
  m_tfinal=0.0;
  if(DEBUG) std::cout<<Pulses.size()<<std::endl;
  std::vector<GRBobsPulse*>::iterator pos;
  for(pos = Pulses.begin(); pos !=  Pulses.end(); ++pos)
    {
      m_tfinal=TMath::Max(m_tfinal,(*pos)->GetEndTime());      
      if(DEBUG) (*pos)->Print();
    }
  
  m_tbin = TMath::Max(10,int(m_tfinal/s_TimeBinWidth));
  m_tbin = TMath::Min(10000,m_tbin);
  s_TimeBinWidth = m_tfinal/m_tbin;
  //  double dt = m_tfinal/(m_tbin-1);
  gDirectory->Delete("Nv");
  m_Nv = new TH2D("Nv","Nv",m_tbin,0.,m_tfinal,Ebin, e);
  
  std::string name;
  GetUniqueName(m_Nv,name);
  m_Nv->SetName(name.c_str());
  
  double t = 0.0;
  for(int ti = 0; ti<m_tbin; ti++)
    {
      t = ti*s_TimeBinWidth;
      for(int ei = 0; ei < Ebin; ei++)
	{
	  double nv = 0.0;//m_Nv->GetBinContent(ti+1, ei+1);
	  
	  for(pos = Pulses.begin(); pos !=  Pulses.end(); ++pos)
	    {	
	      nv += (*pos)->PulseShape(t,e[ei]);
	    }
	  m_Nv->SetBinContent(ti+1, ei+1, nv);
	  // [ph/(cm² s keV)]
	}
    }

  TH2D *nph = Nph(m_Nv); //ph/cm²
  // Scale AT BATSE FLUENCE:
  double norm=0;
  
  if(m_params->GetNormType()==0)
    {
      double BATSEfluence = m_params->GetFluence();
      std::cout<<" Scale at BATSE fluence:" << BATSEfluence <<std::endl;
      
      int ei1 = nph->GetYaxis()->FindBin(BATSE1);
      int ei2 = nph->GetYaxis()->FindBin(BATSE5);
      double F=0.0;
      double en;
      for (int ei = ei1; ei<=ei2; ei++)
	{
	  en   = nph->GetYaxis()->GetBinCenter(ei);
	  for(int ti = 1; ti<=m_tbin; ti++)
	    {
	      F+= nph->GetBinContent(ti, ei) * en;//[keV/cm²]
	    }  
	}
      F*=1.0e-3/(erg2meV); //erg/cm²
      // IMPORTANT m_Nv has to  be in [ph/(m² s keV)]
      norm = 1.0e4 * BATSEfluence/F;
    }
  else
    {
      //////////////////////////////////////////////////
      //SCALE AT BATSE PEAK FLUX:
      double BATSEPeakFlux = m_params->GetPeakFlux();
      std::cout<<" Scale at BATSE PeakFlux:" << BATSEPeakFlux <<std::endl;
      
      int ei1 = nph->GetYaxis()->FindBin(BATSE2);
      int ei2 = nph->GetYaxis()->FindBin(BATSE4);
      double PF=0.0;
      for (int ei = ei1; ei<=ei2; ei++)
	{
	  for(int ti = 1; ti<=m_tbin; ti++)
	    {
	      PF= TMath::Max(PF,nph->GetBinContent(ti, ei)/s_TimeBinWidth);//[ph/cm^2/s]
	    }
	}
      norm = 1.0e4 * BATSEPeakFlux/PF;
    }
  
  m_Nv->Scale(norm);
  Pulses.erase(Pulses.begin(), Pulses.end());
  delete e;
  delete nph;
  return m_Nv;
}
//////////////////////////////////////////////////
TH2D *GRBobsSim::Nph(const TH2D *Nv)
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
      for(int ti = 0; ti<m_tbin; ti++)
	{
	  Nph->SetBinContent(ti+1, ei+1, 
			     Nph->GetBinContent(ti+1, ei+1)*dei*deltat); //[ph/(m²)]  
	}   
    }
  return Nph;
}

//////////////////////////////////////////////////
void GRBobsSim::SaveNv()
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
  sprintf(root_name,"grbobs_%d.root",(int)m_params->GetGRBNumber());
  std::cout<<" Saving "<<root_name<<std::endl;
  TFile mod(root_name,"RECREATE");
  std::string name = m_Nv->GetName();
  m_Nv->SetName("Nv"); // I need a default name.
  m_Nv->Write();
  mod.Close();
  m_Nv->SetName(name.c_str());
  
};

//////////////////////////////////////////////////
void GRBobsSim::SaveGBMDefinition(TString GRBname, double ra, double dec, double theta, double phi, double tstart)
{
  TString name = "GRBOBS_";
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

double LogFBOBS(double *var, double *par)
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
  
  //  if((a-b)>=1.0e-4)
  //    {
  LogH   = log10(a-b) + LogE0;  
  LogC   = (a-b) * (LogH-2.0)-loge*pow(10.0,LogH-LogE0);
  //    }
  //  else
  //    {
  //      a=b+1.0e-4;
  //      LogH  = -4.0 + LogE0;  
  //      LogC  = 1.0e-4 * ((LogH-2.0)- loge);
  //    }
  
  if(LogE <= LogH) 
    return      LogNT + a * (LogE-2.0) - pow(10.0,LogE-LogE0)*loge; 
  return LogC + LogNT + b * (LogE-2.0); // cm^(-2) s^(-1) keV^(-1) 
}



void GRBobsSim::GetGBMFlux(TString GRBname)
{
  //  double *e = new double[Ebin +1];
  //  for(int i = 0; i<=Ebin; i++)
  //    {
  //      e[i] = emin*pow(de,1.0*i); //keV
  //    }
  //////////////////////////////////////////////////
  // GBM Spectrum:
  TF1 band("grb_f",LogFBOBS,log10(emin), 4.0, 4); 
  band.SetParNames("a","b","logE0","Log10(Const)");
 
  band.SetParLimits(0,-2.0 , 2.0); // a
  band.SetParLimits(1,-3.0 , -0.001); // b-a; b < a -> b-a < 0 !
  band.SetParLimits(2,log10(emin),4.0);
  //  band.SetParLimits(3,-3.,3.);

  double a,b,E0,Const;
  TH1D GBM("GBM","GBM",Ebin,log10(emin),log10(emax));
  GBM.SetMinimum(5);
  GBM.SetMaximum(-5);
  
  double t    = 0;
  double dt   = m_Nv->GetXaxis()->GetBinWidth(1);
  double tbin = m_Nv->GetXaxis()->GetNbins();
  
  TString name = "GRBOBS_";
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
	  //	  gPad->SetLogx();
	  //	  gPad->SetLogy();
	  gPad->Update();
	  TString gbmFlux= "GBMFlux";
	  gbmFlux+=ti;
	  gbmFlux+=".gif";
	  if(ti%10==0) gPad->Print(gbmFlux);
	}
    }   
  os.close();
  //////////////////////////////////////////////////
  //  delete[] e;
}
