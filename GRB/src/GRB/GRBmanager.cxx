#include "GRBmanager.h"
#include <iostream>
#include <fstream>
#include "flux/SpectrumFactory.h" 
#include "astro/GPS.h"
#include "astro/SkyDir.h"
#include "astro/PointingTransform.h"
#include "astro/JulianDate.h"
#include "astro/EarthOrbit.h"

ISpectrumFactory &GRBmanagerFactory() 
{
  static SpectrumFactory<GRBmanager> myFactory;
  return myFactory;
}


GRBmanager::GRBmanager(const std::string& params)
  : m_params(params)
{
  m_Nbursts=0;
  paramFile = "$(GRBROOT)/src/test/GRBParam.txt";
  facilities::Util::expandEnvVar(&paramFile);
  
  m_startTime   = TMath::Max(0.,parseParamList(params,0));
  m_timeToWait  = TMath::Max(0.,parseParamList(params,1));
  m_enph        = TMath::Max(0.,parseParamList(params,2))*1000.0; //keV

  m_par = new Parameters();
  //////////////////////////////////////////////////
  GenerateGRB();  
  //////////////////////////////////////////////////
  if(m_enph<=0) m_enph=cst::enph;
}

GRBmanager::~GRBmanager() 
{  
  //  std::cout<<"~GRBmanager() "<<std::endl;
  delete m_par;
  delete m_GRB;
  delete m_spectrum;
}

//return flux, given a time
double GRBmanager::flux(double time) const
{


  double flux;	  
  if(time <= m_startTime || (time > m_endTime)) 
    flux = 0.0;
  else 
    flux = m_spectrum->flux(time-m_startTime,m_enph);
  return flux;
}

TString GRBmanager::GetGRBname(double time)
{
  ////DATE AND GRBNAME
  astro::EarthOrbit m_EarthOrbit;
  using astro::JulianDate;
  JulianDate  JD = m_EarthOrbit.dateFromSeconds(time);
  int An,Me,Gio;
  m_Rest=((int) m_startTime % 86400) + m_startTime - floor(m_startTime);
  m_Frac=(m_Rest/86400.);
  int FracI=(int)(m_Frac*100.0);
  double utc;
  JD.getGregorianDate(An,Me,Gio,utc);

  TString GRBname="";

  An-=2000;
  
  if(An<10) 
    {
      GRBname+="0";
      GRBname+=An;
    }
  else
    {
      GRBname+=An;
    }
  if(Me<10) 
    {
      GRBname+="0";
      GRBname+=Me;
    }
  else
    {
      GRBname+=Me;
    }
  if(Gio<10) 
    {
      GRBname+="0";
      GRBname+=Gio;
    }
  else
    {
      GRBname+=Gio;
    }
  
  if(FracI==100) 
    GRBname+=FracI;
  else if (FracI>10) 
    {
      GRBname+="0";
      GRBname+=FracI;
    }
  else 
    {
      GRBname+="00";
      GRBname+=FracI;
    }
  
  //std::cout<<"GENERATE GRB ("<<GRBname<<")"<<std::endl;
  
  //..................................................  //
  return GRBname;
}


void GRBmanager::GenerateGRB()
{
  TString GRBname = GetGRBname(m_startTime);
  //////////////////////////////////////////////////
  m_Nbursts++;
  m_par->ComputeParametersFromFile(paramFile,m_Nbursts);
  m_GRB      = new  GRBSim(m_par);
  m_spectrum = new  SpectObj(m_GRB->Fireball(),0);
  m_spectrum->SetAreaDetector(EventSource::totalArea());
  //////////////////////////////////////////////////
  m_endTime   = m_startTime + m_GRB->Tmax();
  m_nextBurst = m_endTime   + m_par->rnd->Exp(m_timeToWait);
  m_fluence   = m_GRB->GetFluence();
  m_GRBnumber = m_GRB->GetGRBNumber();

  m_theta = -90.0;
  double FOV=20.0; // degrees above the XY plane.
  while(m_theta<FOV)
    {
      m_par->SetGalDir(-200,-200); //this generates random direction in the sky
      m_GalDir = m_par->GetGalDir();      
      m_l = m_GalDir.first;
      m_b = m_GalDir.second;
      
      astro::SkyDir sky(m_l,m_b,0);
      HepVector3D skydir=sky.dir();
      HepRotation rottoglast = GPS::instance()->transformToGlast(m_startTime,GPS::CELESTIAL);
      HepVector3D scdir = rottoglast * skydir;
      m_ra    = sky.ra();
      m_dec   = sky.dec();
      m_theta = 90. - scdir.theta()*180.0/M_PI; // theta=0 -> XY plane, theta=90 -> Z
      m_phi   = scdir.phi()*180.0/M_PI;
    }

  /////GRB GRNERATION////////////////////////////////
  /*
    TString name = "GRB_";
    name+=GRBname; 
    name+="_PAR.txt";
  */
  
  std::ofstream os;
  if(m_Nbursts==1) 
    os.open("grb_generated.txt",std::ios::out);
  else 
    os.open("grb_generated.txt",std::ios::app);
  
  os<<m_GRBnumber<<" "<<GRBname<<" "<<m_startTime<<" "<<m_endTime<<" "<<m_l<<" "<<m_b<<" "<<m_theta<<" "<<m_phi<<" "<<m_fluence<<" "<<std::endl;
  os.close();
  std::cout<<"Physical Model GRB"<<GRBname<<" t start "<<m_startTime<<", tend "<<m_endTime
	   <<" l,b = "<<m_l<<", "<<m_b<<" elevation,phi(deg) = "<<m_theta<<", "<<m_phi<<" Fluence = "<<m_fluence<<std::endl;
  
  if(m_par->GenerateGBM())
    {
      m_GRB->SaveGBMDefinition(GRBname,m_ra,m_dec,m_theta,m_phi,m_startTime);
      m_GRB->GetGBMFlux(GRBname);
    }
  m_startTime = m_nextBurst;
}

double GRBmanager::interval(double time)
{  
  double inte;  
  
  if(time <= m_startTime) 
    inte = m_startTime - time + m_spectrum->interval(0.0,m_enph);
  else if (time<m_endTime)
    inte = m_spectrum->interval(time - m_startTime,m_enph);
  else 
    {
      delete m_GRB;
      delete m_spectrum;
      GenerateGRB(); 
      inte = m_startTime-time + m_spectrum->interval(0.0,m_enph);
    }
  inte = TMath::Min(inte,m_nextBurst-time);
  //  std::cout<<"GRBmanager interval : "<<inte<<std::endl;
  return inte;
}

double GRBmanager::energy(double time)
{
  double energy=m_spectrum->energy(time-m_startTime,m_enph)*1.0e-3; //MeV
  //  std::cout<<"GRBmanager energy "<<energy<<std::endl;
  return energy;
}

double GRBmanager::parseParamList(std::string input, int index)
{
  std::vector<double> output;
  unsigned int i=0;
  for(;!input.empty() && i!=std::string::npos;){
    double f = ::atof( input.c_str() );
    output.push_back(f);
    i=input.find_first_of(",");
    input= input.substr(i+1);
  } 
  if(index>=(int) output.size()) return 0.0;
  return output[index];
}


