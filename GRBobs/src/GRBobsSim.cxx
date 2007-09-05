#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm> //sort
#include <string>

#include <stdexcept>

#include "GRBobs/GRBobsConstants.h"
#include "GRBobs/GRBobsengine.h"
#include "GRBobs/GRBobsSim.h"
#include "GRBobs/GRBobsPulse.h"

#include "TFile.h"
#include "TCanvas.h"
#define DEBUG 0
#define APPLY_REDSHIFT 0

using namespace ObsCst;
//////////////////////////////////////////////////
double *GRBobsSim::ComputeEnergyBins(int &Nbins)
{
  NaIEnergyGrid_Vector.erase(NaIEnergyGrid_Vector.begin(),NaIEnergyGrid_Vector.end());
  BGOEnergyGrid_Vector.erase(BGOEnergyGrid_Vector.begin(),BGOEnergyGrid_Vector.end());

  std::vector<double> EnergyGrid;
  std::string path = ::getenv("GRBOBSROOT");
  std::string NaIpath(path+"/GRBobs/NaI_energy_grid.dat");
  std::string BGOpath(path+"/GRBobs/BGO_energy_grid.dat");
  //  std::cout<<NaIpath<<" "<<BGOpath<<std::endl;
  std::ifstream NaIEnergyGrid(NaIpath.c_str());

  if(!NaIEnergyGrid.is_open())
    {
      std::cout<<"Unable to open "+NaIpath<<std::endl;
      throw std::runtime_error("Unable to open "+NaIpath);
    }
  std::ifstream BGOEnergyGrid(BGOpath.c_str());
  if(!BGOEnergyGrid.is_open())
    {
      std::cout<<"Unable to open "+NaIpath<<std::endl;
      throw std::runtime_error("Unable to open "+BGOpath);
    }
  
  int NaI_number_of_energy_bins;
  int BGO_number_of_energy_bins;

  double energy_low;
  
  NaIEnergyGrid>>NaI_number_of_energy_bins;

  
  for(int i = 0; i<NaI_number_of_energy_bins; i++)
    {
      NaIEnergyGrid >> energy_low;
      NaIEnergyGrid_Vector.push_back(energy_low);
      EnergyGrid.push_back(energy_low);
    }
  
  
  BGOEnergyGrid>>BGO_number_of_energy_bins;
  
  for(int i = 0; i<BGO_number_of_energy_bins; i++)
    {
      BGOEnergyGrid >> energy_low;
      EnergyGrid.push_back(energy_low);
      BGOEnergyGrid_Vector.push_back(energy_low);
    }
  // LAT_number_of_energy_bins is in ObsCst!!!
  double delta_e = pow(emax/energy_low,1./LAT_number_of_energy_bins);
  int energy_bin = 1;
  while(energy_bin < LAT_number_of_energy_bins)
    {
      EnergyGrid.push_back(energy_low*pow(delta_e,energy_bin++));
    }
  EnergyGrid.push_back(emax);
  //..................................................//
  
  std::sort(EnergyGrid.begin(),EnergyGrid.end());
  
  /*
    for(std::vector<double>::iterator pos=EnergyGrid.begin(); pos!=EnergyGrid.end();pos++)
    std::cout<<(*pos)<<std::endl;
  */
  std::vector<double>::iterator pos=EnergyGrid.begin(); 
  while (pos!=EnergyGrid.end())
    {
      if(*(pos)==*(pos+1)) EnergyGrid.erase(pos);
      else
	pos++;
    }
  //////////////////////////////////////////////////
  Nbins = EnergyGrid.size();
  double *e = new double[Nbins];
  Nbins--;
  for(int i=0; i<= Nbins; i++)
    {
      e[i]= EnergyGrid[i];
    }
  return e;
}


//////////////////////////////////////////////////

GRBobsSim::GRBobsSim(GRBobsParameters *params)
  : m_params(params)
{
  m_GRBengine = new GRBobsengine(params);
}

void GRBobsSim::GetUniqueName(void *ptr, std::string & name)
{
  std::ostringstream my_name;
  my_name << reinterpret_cast<int> (ptr);
  name = my_name.str();
  gDirectory->Delete(name.c_str());
  reinterpret_cast<TH1*> (ptr)->SetDirectory(0);
}

TH2D* GRBobsSim::MakeGRB()
{
  double z = m_params->GetRedshift();
  double duration =  m_params->GetDuration();
  int Ebin;
  double *energies = ComputeEnergyBins(Ebin);
  double s_TimeBinWidth=TimeBinWidth;
  std::vector<GRBobsPulse*> Pulses = m_GRBengine->CreatePulsesVector();
  m_tfinal=0.0;
  if(DEBUG) std::cout<<Pulses.size()<<std::endl;
  std::vector<GRBobsPulse*>::iterator pos;
  
  //  for(pos = Pulses.begin(); pos !=  Pulses.end(); ++pos)
  //    {
  //      m_tfinal= TMath::Max(m_tfinal,(*pos)->GetEndTime());      
  //      if(DEBUG) (*pos)->Print();
  //    }
  
  if(APPLY_REDSHIFT) m_tfinal=(1.0+z) * duration; // duration is in the intrinsic frame;
  else m_tfinal= duration; // duration is in the observer frame;
  
  m_tbin = TMath::Max(10,int(m_tfinal/s_TimeBinWidth));
  //  m_tbin = TMath::Min(10000,m_tbin);
  s_TimeBinWidth = m_tfinal/m_tbin;
  
  gDirectory->Delete("Nv");
  m_Nv = new TH2D("Nv","Nv",m_tbin,0.,m_tfinal,Ebin, energies);
  
  std::string name;
  GetUniqueName(m_Nv,name);
  m_Nv->SetName(name.c_str());

  double Fssc_Fsyn = m_params->GetFssc_Fsyn();
  double Essc_Esyn = m_params->GetEssc_Esyn();
  

  for(int ti = 0; ti<m_tbin; ti++)
    {
      double t = m_Nv->GetXaxis()->GetBinCenter(ti+1);
      for(int ei = 0; ei < Ebin; ei++)
	{
	  double e  = m_Nv->GetYaxis()->GetBinCenter(ei+1);
	  double nv = 0.0;
	  
	  for(pos = Pulses.begin(); pos !=  Pulses.end(); ++pos)
	    {
	      if(APPLY_REDSHIFT) nv += (*pos)->PulseShape(t/(1.+z),e*(1.+z)); //t/(1.+z) and e*(1.+z) are intrinsic
	      else nv += (*pos)->PulseShape(t,e); //t and e are observed
	      if(Essc_Esyn*Fssc_Fsyn>0.0) // Add the ssc component
		{
		  if(APPLY_REDSHIFT) nv += Fssc_Fsyn*(*pos)->PulseShape(t/(1.+z),e*(1.+z)/Essc_Esyn)/(Essc_Esyn*Essc_Esyn); //t/(1.+z) and e*(1.+z) are intrinsic
		  else nv += Fssc_Fsyn*(*pos)->PulseShape(t,e/Essc_Esyn)/(Essc_Esyn*Essc_Esyn); //t and e are observed
		}
	    }
	  m_Nv->SetBinContent(ti+1, ei+1, nv);
	  // [ph/(cm² s keV)]
	}
    }
  //////////////////////////////////////////////////
  TH2D *nph = Nph(m_Nv); //ph/cm²
  //////////////////////////////////////////////////
  // Scale AT BATSE FLUENCE:
  double norm=0;
  
  if(m_params->GetNormType()=='F')
    {
      double BATSEfluence = m_params->GetFluence();
      if(DEBUG) std::cout<<" Scale at BATSE fluence!:" << BATSEfluence <<std::endl;
      
      int ei1 = nph->GetYaxis()->FindBin(BATSE2);
      int ei2 = nph->GetYaxis()->FindBin(BATSE4);
      double F=0.0;
      for (int ei = ei1; ei<=ei2; ei++)
	{
	  double e   = nph->GetYaxis()->GetBinCenter(ei);
	  for(int ti = 1; ti<=m_tbin; ti++)
	    {
	      F+= nph->GetBinContent(ti, ei) * e;//[keV/cm²]
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
      if(DEBUG) std::cout<<" Scale at BATSE PeakFlux!:" << BATSEPeakFlux <<std::endl;
      
      int ei1 = nph->GetYaxis()->FindBin(BATSE2);
      int ei2 = nph->GetYaxis()->FindBin(BATSE4);
      
      //////////////////////////////////////////////////
      double AccumulationTime = 0.256;//s_TimeBinWidth; (BONNELL et al. 1997, Apj.)
      double PF=0.0;
      double acc_time=0.0;
      double PF_acc=0.0;
      for(int ti = 1; ti<=m_tbin; ti++)
	{
	  double PF_t=0.0;
	  acc_time += s_TimeBinWidth;
	  for (int ei = ei1; ei<=ei2; ei++)
	    {
	      PF_t += nph->GetBinContent(ti, ei); //[ph]
	    }
	  
	  if(acc_time <= AccumulationTime && ti < m_tbin)
	    {
	      PF_acc += PF_t;
	    }
	  else
	    {
	      PF = TMath::Max(PF,PF_acc/acc_time); //[ph/s]
	      acc_time = 0.0;
	      PF_acc   = 0.0;
	    }
	  //PF = TMath::Max(PF,PF_t); //[ph]
	}
      
      norm = 1.0e4 * BATSEPeakFlux/PF;
    }
  
  m_Nv->Scale(norm);
  Pulses.erase(Pulses.begin(), Pulses.end());
  for(int i =0; i<(int) Pulses.size();i++) delete Pulses[i];
  delete[] energies;
  delete nph;
  return m_Nv;
}

//////////////////////////////////////////////////
TH2D* GRBobsSim::MakeGRB_ExtraComponent(double duration, double LATphotons)
{
  int Ebin;
  double *energies = ComputeEnergyBins(Ebin);
  const int tbin = 1000;
  const double dt = duration/(1.0*tbin);
  
  gDirectory->Delete("NvEC");
  m_NvEC = new TH2D("NvEC","NvEC",tbin,0.,duration,Ebin, energies);
  
  std::string name;
  GetUniqueName(m_NvEC,name);
  m_NvEC->SetName(name.c_str());
  
  for(int ti = 0; ti<tbin; ti++)
    {
      double t  = ti*dt;
      double I0 = pow(duration/(9.*t+duration),2.0);
      for(int ei = 0; ei < Ebin; ei++)
	{
	  double e  = m_Nv->GetYaxis()->GetBinCenter(ei+1);
	  double I = I0*pow(e,-1.0);
	  m_NvEC->SetBinContent(ti+1,ei+1,I); // [ph/(cm² s keV)]
	}
    }
  
  int ei1 = m_NvEC->GetYaxis()->FindBin(LAT1);
  int ei2 = m_NvEC->GetYaxis()->FindBin(LAT2);
  double norm=0.0;  
  for(int ti = 1; ti<=tbin; ti++)
    {
      for (int ei = ei1; ei<=ei2; ei++)
	{
	  double de = m_NvEC->GetYaxis()->GetBinWidth(ei);
	  norm += m_NvEC->GetBinContent(ti, ei) * dt * de; // ph/(cm²)
	}
    }
  //  std::cout<<LATphotons<<" "<<norm<<std::endl;
  m_NvEC->Scale(LATphotons/norm);
  delete[] energies;
  return m_NvEC;
}
//////////////////////////////////////////////////

TH2D* GRBobsSim::CutOff(TH2D *Nv, double E_CO)
{
  
  if (E_CO==0) return Nv;
  int tbin = Nv->GetXaxis()->GetNbins();
  int Ebin = Nv->GetYaxis()->GetNbins();
  for(int ti = 0; ti<tbin; ti++)
    {
      for(int ei = 0; ei < Ebin; ei++)
	{
	  double nv = Nv->GetBinContent(ti+1,ei+1); // [ph/(cm² s keV)]
	  double e = Nv->GetYaxis()->GetBinCenter(ei+1)*1e-6; // (GeV)
	  double suppression = exp(-e/E_CO);
	  Nv->SetBinContent(ti+1,ei+1,suppression*nv); // [ph/(cm² s keV)]
	}
    }
  return Nv;
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
  int Ebin      = Nv->GetYaxis()->GetNbins();
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
void GRBobsSim::SaveNvEC()
{
  
  m_NvEC->SetXTitle("Time [s]");
  m_NvEC->SetYTitle("Energy [keV]");
  m_NvEC->SetZTitle("N_{v} [ph/m^2/s/keV]");
  m_NvEC->GetXaxis()->SetTitleOffset(1.5);
  m_NvEC->GetYaxis()->SetTitleOffset(1.5);
  m_NvEC->GetZaxis()->SetTitleOffset(1.2);
  m_NvEC->GetXaxis()->CenterTitle();
  m_NvEC->GetYaxis()->CenterTitle();
  m_NvEC->GetZaxis()->CenterTitle();
  
  char root_name[100];
  sprintf(root_name,"grbobs_%d_EC.root",(int)m_params->GetGRBNumber());
  if(DEBUG) std::cout<<" Saving "<<root_name<<std::endl;
  TFile mod(root_name,"RECREATE");
  std::string name = m_NvEC->GetName();
  m_NvEC->SetName("Nv"); // I need a default name.
  m_NvEC->Write();
  mod.Close();
  m_NvEC->SetName(name.c_str());
};

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
  if(DEBUG) std::cout<<" Saving "<<root_name<<std::endl;
  TFile mod(root_name,"RECREATE");
  std::string name = m_Nv->GetName();
  m_Nv->SetName("Nv"); // I need a default name.
  m_Nv->Write();
  mod.Close();
  m_Nv->SetName(name.c_str());

};

//////////////////////////////////////////////////
void GRBobsSim::SaveGBMDefinition(std::string GRBname, double ra, double dec, double theta, double phi, double tstart)
{
  std::string name = "GRBOBS_";
  name+=GRBname; 
  name+=".DEF";
  std::ofstream os(name.c_str(),std::ios::out);
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

/*
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
*/


void GRBobsSim::GetGBMFlux(std::string GRBname)
{
  //  m_Nv has to  be in [ph/(m² s keV)]
  if(DEBUG) std::cout<<" NaI channels: "<<NaIEnergyGrid_Vector.size()<<" BGO channels: "<<BGOEnergyGrid_Vector.size()<<std::endl;
  //  double t    = 0;
  //  double dt   = m_Nv->GetXaxis()->GetBinWidth(1);
  double tbin = m_Nv->GetXaxis()->GetNbins();
  //  double ebin = m_Nv->GetYaxis()->GetNbins();
  
  std::string name_NaI = "NaI_GRBOBS_";
  name_NaI += GRBname; 
  name_NaI += ".lc";

  std::string name_BGO = "BGO_GRBOBS_";
  name_BGO += GRBname; 
  name_BGO += ".lc";
  
  std::ofstream file_NaI(name_NaI.c_str(),std::ios::out);
  std::ofstream file_BGO(name_BGO.c_str(),std::ios::out);

  file_NaI<<"TimeBins= ";
  file_BGO<<"TimeBins= ";

  file_NaI<<tbin<<std::endl;
  file_BGO<<tbin<<std::endl;

  TH1D NaISpectrum("NaISpectrum","NaISpectrum",NaIEnergyGrid_Vector.size()-1,&NaIEnergyGrid_Vector[0]); 
  TH1D BGOSpectrum("BGOSpectrum","BGOSpectrum",BGOEnergyGrid_Vector.size()-1,&BGOEnergyGrid_Vector[0]); 
  NaISpectrum.SetLineColor(2);
  BGOSpectrum.SetLineColor(4);

#ifndef WIN32 // THB: Avoid need to link TCanvas on windows 
  if(DEBUG)
    {
      TCanvas *GBMCanvas;
      GBMCanvas = new TCanvas("GBMCanvas","GBMCanvas",500,400);
      GBMCanvas->SetLogx();
      GBMCanvas->SetLogy();
    }
#endif
  
  
  
  for(int ti=1; ti<=tbin; ti++)
    {
      for(unsigned int ei = 0; ei<NaIEnergyGrid_Vector.size()-1; ei++)
	{
	  double sp=0.0;
      	  double e_low  = NaIEnergyGrid_Vector[ei];
	  double e_high = NaIEnergyGrid_Vector[ei+1];
	  
	  int    ei_low = m_Nv->GetYaxis()->FindBin(e_low);
	  int    ei_high = m_Nv->GetYaxis()->FindBin(e_high);
	  
	  int diff = ei_high - ei_low;
	  for(int i=0;i<diff;i++)
	    {
	      sp+=(m_Nv->GetBinContent(ti,ei_low+i)*m_Nv->GetYaxis()->GetBinWidth(ei_low+i)); //ph/m^2/s 
	    }
	  NaISpectrum.SetBinContent(ei,sp/(e_high-e_low)*1e-4); //[ph/(cm² s keV)]
	}
      //////////////////////////////////////////////////
      for(unsigned int ei = 0; ei<BGOEnergyGrid_Vector.size()-1; ei++)
	{
	  double sp=0.0;
	  
      	  double e_low  = BGOEnergyGrid_Vector[ei];
	  double e_high = BGOEnergyGrid_Vector[ei+1];
	  int    ei_low = m_Nv->GetYaxis()->FindBin(e_low);
	  int    ei_high = m_Nv->GetYaxis()->FindBin(e_high);
	  
	  int diff = ei_high - ei_low;
	  for(int i=0;i<diff;i++)
	    {
	      sp+=(m_Nv->GetBinContent(ti,ei_low+i)*m_Nv->GetYaxis()->GetBinWidth(ei_low+i)); //ph/m^2/s ;
	    }
	  BGOSpectrum.SetBinContent(ei,sp/(e_high-e_low)*1e-4); //[ph/(cm² s keV)]
	}
      
#ifndef WIN32 // THB: Avoid need to link TCanvas on windows 
      if(DEBUG)
	{
	  NaISpectrum.Draw();
	  BGOSpectrum.Draw("same");
	  gPad->Update();
	}
#endif
      for(unsigned int ei = 1; ei<NaIEnergyGrid_Vector.size(); ei++) 
	file_NaI << std::setprecision(5) << NaISpectrum.GetBinContent(ei)<<"\t";
      file_NaI <<"\n";
      for(unsigned int ei = 1; ei<BGOEnergyGrid_Vector.size(); ei++) 
	file_BGO << std::setprecision(5) << BGOSpectrum.GetBinContent(ei)<<"\t";
      file_BGO <<"\n";
    }
  file_NaI.close();
  file_BGO.close();
  
}
