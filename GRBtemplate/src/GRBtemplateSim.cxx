#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <stdexcept>
#include <cmath>
#include "GRBtemplate/GRBtemplateSim.h"
#include "facilities/commonUtilities.h"

#include "TFile.h"
#include "TCanvas.h"
#define DEBUG 0

using namespace TmpCst;

TH1D *projectHistogram(TH1D *Nv, std::vector<double> energyBins)
{
  int N2 = energyBins.size();
  gDirectory->Delete("GBM");
  TH1D *GBM = new TH1D("GBM","GBM",energyBins.size()-1,&energyBins[0]);

  int N1 = Nv->GetNbinsX();

  for (int i=0;i<N2;i++)
    {
      double e1 = energyBins[i];
      double e2 = energyBins[i+1];
      double x   = (e2+e1)/2.;
      
      int    I0      = Nv->FindBin(e1);
      int    I1      = Nv->FindBin(e2);

      double X,X1,X2,Y1,Y2;
      double y=0.0;      
      if(N2>=N1)
	{
	  X     =  Nv->GetBinCenter(I0);
	  if(x<X)
	    {
	      if(I0==1)
		{
		  X1     =  Nv->GetBinCenter(I0);
		  X2     =  Nv->GetBinCenter(I0+1);
		  Y1     =  Nv->GetBinContent(I0);
		  Y2     =  Nv->GetBinContent(I0+1);
		}
	      else
		{
		  X1     =  Nv->GetBinCenter(I0-1);
		  X2     =  Nv->GetBinCenter(I0);
		  Y1     =  Nv->GetBinContent(I0-1);
		  Y2     =  Nv->GetBinContent(I0);
		}
	    }
	  else
	    {
	      if(I0==N1)
		{
		  X1     =  Nv->GetBinCenter(I0-1);
		  X2     =  Nv->GetBinCenter(I0);
		  Y1     =  Nv->GetBinContent(I0-1);
		  Y2     =  Nv->GetBinContent(I0);
		}
	      else
		{
		  X1     =  Nv->GetBinCenter(I0);
		  X2     =  Nv->GetBinCenter(I0+1);
		  Y1     =  Nv->GetBinContent(I0);
		  Y2     =  Nv->GetBinContent(I0+1);
		}
	    }
	  y = Y1+(Y2-Y1)/(X2-X1)*(x-X1);
	}
      else
	{
	  int diff = I1-I0;
	  y=0.0;
	  if(diff==0)
	    {
	      double de12 = e2 - e1;
	      y+= Nv->GetBinContent(I0)*de12;
	    }
	  else if(diff>0)
	    {
	      double de1 = Nv->GetBinLowEdge(I0)+ Nv->GetBinWidth(I0) - e1;
	      y+= Nv->GetBinContent(I0)*de1;
	      double de2 = e2 - Nv->GetBinLowEdge(I1);
	      y+= Nv->GetBinContent(I1)*de2;
	      if(diff>1)
		for(int j=1;j<diff;j++)
		  {
		    y+=(Nv->GetBinContent(I0+j)*Nv->GetBinWidth(I0+j)); //ph/m^2/s 
		  }
	    }
	}
      y/=(e2-e1); //[ph/(cm² s keV)]
      GBM->SetBinContent(i+1,y);
    }
  return GBM;
}



void GRBtemplateSim::ComputeEnergyBins()
{
  NaIEnergyGrid_Vector.erase(NaIEnergyGrid_Vector.begin(),NaIEnergyGrid_Vector.end());
  BGOEnergyGrid_Vector.erase(BGOEnergyGrid_Vector.begin(),BGOEnergyGrid_Vector.end());
  
  std::string NaIpath(facilities::commonUtilities::joinPath(facilities::commonUtilities::getDataPath("GRBtemplate"),"NaI_energy_grid.dat"));
  std::string BGOpath(facilities::commonUtilities::joinPath(facilities::commonUtilities::getDataPath("GRBtemplate"),"BGO_energy_grid.dat"));
  std::cout<<NaIpath<<" *** "<<BGOpath<<std::endl;
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
    }
  
  
  BGOEnergyGrid>>BGO_number_of_energy_bins;
  
  for(int i = 0; i<BGO_number_of_energy_bins; i++)
    {
      BGOEnergyGrid >> energy_low;
      BGOEnergyGrid_Vector.push_back(energy_low);
    }
}

//////////////////////////////////////////////////
GRBtemplateSim::GRBtemplateSim(std::string InputFileName)
  : m_InputFileName(InputFileName)
{
  
}

void GRBtemplateSim::GetUniqueName(void *ptr, std::string & name)
{
  std::ostringstream my_name;
  my_name << reinterpret_cast<int> (ptr);
  name = my_name.str();
  gDirectory->Delete(name.c_str());
  reinterpret_cast<TH1*> (ptr)->SetDirectory(0);
}

TH2D* GRBtemplateSim::MakeGRB()
{
  ifstream iFile(m_InputFileName.c_str());
  char dummy[100];
  iFile>>dummy>>m_EnergyBins;
  iFile>>dummy>>m_emin;
  iFile>>dummy>>m_emax;
  iFile>>dummy>>m_TimeBins;
  iFile>>dummy>>m_TimeBinWidth;

  std::cout<<"Template from: "<<m_InputFileName<<std::endl;
  std::cout<<"EnergyBins= "<<m_EnergyBins<<std::endl;
  std::cout<<"MinimumEnergy= "<<m_emin<<std::endl;
  std::cout<<"MaximumEnergy= "<<m_emax<<std::endl;
  std::cout<<"TimeBins= "<<m_TimeBins<<std::endl;
  std::cout<<"TimeBinWidth= "<<m_TimeBinWidth<<std::endl;

  double de   = pow(m_emax/m_emin,1.0/m_EnergyBins);  


  double *e  = new double[m_EnergyBins +1];

  for(int i = 0; i<=m_EnergyBins; i++)
    {
     e[i] = m_emin*pow(de,1.0*i); //keV
    }
  //////////////////////////////////////////////////  
  m_tfinal=m_TimeBinWidth * (m_TimeBins-1);
  //////////////////////////////////////////////////
  gDirectory->Delete("Nv");  

  m_Nv = new TH2D("Nv","Nv",m_TimeBins,0.,m_tfinal,m_EnergyBins, e);
  std::string name;
  GetUniqueName(m_Nv,name);
  m_Nv->SetName(name.c_str());
  

  //  double t = 0.0;  
  for(int ti = 0; ti<m_TimeBins; ti++)
    {
      //      t = ti * m_TimeBinWidth;
      for(int ei = 0; ei < m_EnergyBins; ei++)
	{
	  double nv;
	  iFile>>nv; // [ph/(cm² s keV)]
	  m_Nv->SetBinContent(ti+1, ei+1, nv*1e4);// [ph/(m² s keV)]
	}
    }

  TH2D *nph = Nph(m_Nv); //ph/m²
  delete[] e;
  delete nph;
  return m_Nv;
}
//////////////////////////////////////////////////
TH2D *GRBtemplateSim::Nph(const TH2D *Nv)
{
  
  TH2D *Nph = (TH2D*) Nv->Clone(); // [ph/(m² s keV)]  
  std::string name;
  GetUniqueName(Nph,name);
  Nph->SetName(name.c_str());
  
  double dei;
  double deltat = Nv->GetXaxis()->GetBinWidth(1);
  
  for (int ei = 0; ei<m_EnergyBins; ei++)
    {
      dei   = Nv->GetYaxis()->GetBinWidth(ei+1);
      for(int ti = 0; ti<m_TimeBins; ti++)
	{
	  Nph->SetBinContent(ti+1, ei+1, 
			     Nph->GetBinContent(ti+1, ei+1)*dei*deltat); //[ph/(m²)]  
	}   
    }
  return Nph;
}

//////////////////////////////////////////////////
void GRBtemplateSim::SaveNv()
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
  sprintf(root_name,"grbtmp.root");//,(int)m_params->GetGRBNumber());
  std::cout<<" Saving "<<root_name<<std::endl;
  TFile mod(root_name,"RECREATE");
  std::string name = m_Nv->GetName();
  m_Nv->SetName("Nv"); // I need a default name.
  m_Nv->Write();
  mod.Close();
  m_Nv->SetName(name.c_str());
  
};

//////////////////////////////////////////////////
void GRBtemplateSim::SaveGBMDefinition(std::string GRBname, double ra, double dec, double theta, double phi, double tstart)
{
  std::string name = "GRBTMP_";
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

double LogFBTMP(double *var, double *par)
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


void GRBtemplateSim::GetGBMFlux(std::string GRBname)
{
  ComputeEnergyBins();
  //  m_Nv has to  be in [ph/(m² s keV)]
  if(DEBUG) std::cout<<" NaI channels: "<<NaIEnergyGrid_Vector.size()<<" BGO channels: "<<BGOEnergyGrid_Vector.size()<<std::endl;
  double tbin = m_Nv->GetXaxis()->GetNbins();
  
  std::string name_NaI = "NaI_GRBTMP_";
  name_NaI += GRBname; 
  name_NaI += ".lc";
  
  std::string name_BGO = "BGO_GRBTMP_";
  name_BGO += GRBname; 
  name_BGO += ".lc";
  
  std::ofstream file_NaI(name_NaI.c_str(),std::ios::out);
  std::ofstream file_BGO(name_BGO.c_str(),std::ios::out);

  file_NaI<<"TimeBins= ";
  file_BGO<<"TimeBins= ";
  
  file_NaI<<tbin<<std::endl;
  file_BGO<<tbin<<std::endl;
  double *e  = new double[m_EnergyBins +1];
  for (int ei = 0; ei<m_EnergyBins; ei++) 
    e[ei]=m_Nv->GetYaxis()->GetBinLowEdge(ei+1);
  e[m_EnergyBins]=m_Nv->GetYaxis()->GetBinLowEdge(m_EnergyBins)+m_Nv->GetYaxis()->GetBinWidth(m_EnergyBins);
  
  TH1D *Nve =  new TH1D("Nve","Nve",m_EnergyBins, e);
#ifndef WIN32 // THB: Avoid need to link TCanvas on windows 
  if(DEBUG)
    {
      TCanvas *GBMCanvas;
      GBMCanvas = new TCanvas("GBMCanvas","GBMCanvas",1000,800);
      GBMCanvas->SetLogx();
      GBMCanvas->SetLogy();
    }
#endif
  
  for(int ti=1; ti<=tbin; ti++)
    {
      for (int ei = 1; ei<=m_EnergyBins; ei++) Nve->SetBinContent(ei,m_Nv->GetBinContent(ti,ei));
      TH1D *NaISpectrum = projectHistogram(Nve,NaIEnergyGrid_Vector);
      NaISpectrum->SetNameTitle("NaI","NaI");
      TH1D *BGOSpectrum = projectHistogram(Nve,BGOEnergyGrid_Vector);
      NaISpectrum->SetNameTitle("BGO","BGO");
      

      NaISpectrum->SetLineColor(2);
      BGOSpectrum->SetLineColor(4);
#ifndef WIN32 // THB: Avoid need to link TCanvas on windows 
      if(DEBUG)
	{
	  Nve->Draw();
	  NaISpectrum->Draw("same");
	  BGOSpectrum->Draw("same");
	  gPad->Update();
	}
#endif
      for(unsigned int ei = 1; ei<NaIEnergyGrid_Vector.size(); ei++) 
	file_NaI << std::setprecision(5) << NaISpectrum->GetBinContent(ei)<<"\t";
      file_NaI <<"\n";
      for(unsigned int ei = 1; ei<BGOEnergyGrid_Vector.size(); ei++) 
	file_BGO << std::setprecision(5) << BGOSpectrum->GetBinContent(ei)<<"\t";
      file_BGO <<"\n";
      gDirectory->Delete("NaI");
      gDirectory->Delete("BGO");
    }
  file_NaI.close();
  file_BGO.close();
}
