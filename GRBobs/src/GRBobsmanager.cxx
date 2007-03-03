#include "GRBobs/GRBobsmanager.h"
#include <iostream>
#include <ostream>
#include <fstream>
#include "flux/SpectrumFactory.h" 
#include "astro/GPS.h"
#include "astro/SkyDir.h"
#include "astro/PointingTransform.h"
#include "astro/JulianDate.h"
#include "astro/EarthOrbit.h"


#define DEBUG 0 

ISpectrumFactory &GRBobsmanagerFactory() 
{
  static SpectrumFactory<GRBobsmanager> myFactory;
  return myFactory;
}


GRBobsmanager::GRBobsmanager(const std::string& params)
  : m_params(params)
{
  using astro::GPS;
  m_GenerateGBMOutputs = false;

  facilities::Util::expandEnvVar(&paramFile);  
  
  m_startTime       = parseParamList(params,0);
  double duration   = parseParamList(params,1);
  m_fluence         = parseParamList(params,2);
  m_alpha           = parseParamList(params,3);
  m_beta            = parseParamList(params,4);
  m_MinPhotonEnergy = parseParamList(params,5)*1.0e3; //MeV

  if(parseParamList(params,6)!=0) m_GenerateGBMOutputs = true;
  
    
  m_par = new GRBobsParameters();

  m_GRBnumber = (long) floor(65540+m_startTime);

  m_par->SetGRBNumber(m_GRBnumber);
  m_par->SetFluence(m_fluence);
  m_par->SetPeakFlux(m_fluence);
  m_par->SetDuration(duration);
  m_par->SetAlphaBeta(m_alpha,m_beta);
  m_par->SetMinPhotonEnergy(m_MinPhotonEnergy); //keV
  
  m_theta = -90.0;
  double FOV=20.0; // degrees above the XY plane.
  while(m_theta<FOV)
    {
      m_par->SetGalDir(-200,-200); //this generates random direction in the sky
      m_GalDir = m_par->GetGalDir();      
      m_l = m_GalDir.first;
      m_b = m_GalDir.second;
      
      astro::SkyDir sky(m_l,m_b,astro::SkyDir::GALACTIC);
      HepVector3D skydir=sky.dir();
      HepRotation rottoglast = GPS::instance()->transformToGlast(m_startTime,GPS::CELESTIAL);
      HepVector3D scdir = rottoglast * skydir;
      m_ra    = sky.ra();
      m_dec   = sky.dec();
      m_theta = 90. - scdir.theta()*180.0/M_PI; // theta=0 -> XY plane, theta=90 -> Z
      m_phi   = scdir.phi()*180.0/M_PI;
    }

  m_par->SetGRBNumber(m_GRBnumber);
  m_grbGenerated    = false;
  m_grbdeleted      = false;
  //////////////////////////////////////////////////
}

TString GRBobsmanager::GetGRBname()
{
  ////DATE AND GRBNAME
  /// Time from 2005 ,7,18,0.0:
  //  astro::EarthOrbit m_EarthOrbit;
  //  astro::JulianDate JD = m_EarthOrbit.dateFromSeconds(time);
  //  else: 
  astro::JulianDate JD(2007, 1, 1, 0.0); //DC2
  astro::EarthOrbit m_EarthOrbit(JD);
  //  astro::JulianDate JD = m_EarthOrbit.dateFromSeconds(static_cast<double>(mytime));
  //  JD = m_EarthOrbit.mjdFromSeconds(m_startTime);
  JD = m_EarthOrbit.dateFromSeconds(m_startTime);
  int An,Me,Gio;
  m_Rest=((int) m_startTime % 86400) + m_startTime - floor(m_startTime);
  m_Frac=(m_Rest/86400.);
  int FracI=(int)(m_Frac*1000.0);
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
  
  
  if (FracI<10)
    {
      GRBname+="00";
      GRBname+=FracI;
    }
  else if(FracI<100) 
    {
      GRBname+="0";
      GRBname+=FracI;
    }
  else 
    GRBname+=FracI;
  
  if(DEBUG) std::cout<<"GENERATE GRB ("<<GRBname<<")"<<std::endl;
  
  //..................................................  //
  return GRBname;
}

GRBobsmanager::~GRBobsmanager() 
{  
  DeleteGRB();
}

void GRBobsmanager::GenerateGRB()
{
  /////GRB GRNERATION////////////////////////////////
  m_GRB      = new  GRBobsSim(m_par);
  m_spectrum = new  SpectObj(m_GRB->MakeGRB(),0);
  m_spectrum->SetAreaDetector(EventSource::totalArea());
  //////////////////////////////////////////////////
  m_endTime   = m_startTime + m_GRB->Tmax();
  TString GRBname = GetGRBname();

  if(m_GenerateGBMOutputs)
    {
      m_GRB->SaveGBMDefinition(GRBname,m_ra,m_dec,m_theta,m_phi,m_Rest);
      m_GRB->GetGBMFlux(GRBname);
    }
  
  TString name = "GRBOBS_";
  name+=GRBname; 
  name+="_PAR.txt";
  
  std::ofstream os(name,std::ios::out);
  os<<m_startTime<<" "<<m_endTime<<" "<<m_l<<" "<<m_b<<" "<<m_theta<<" "<<m_phi<<" "<<m_fluence<<" "<<" "<<m_alpha<<" "<<m_beta<<std::endl;
  os.close();
  std::cout<<"Phenomen. Model GRB"<<GRBname<<" t start "<<m_startTime<<", tend "<<m_endTime
	   <<" l,b = "<<m_l<<", "<<m_b<<" elevation,phi(deg) = "<<m_theta<<", "<<m_phi<<" Fluence = "<<m_fluence<<std::endl;
  m_grbGenerated=true;
}

void GRBobsmanager::DeleteGRB()
{
  delete m_par;
  delete m_GRB;
  delete m_spectrum;
  m_grbdeleted=true;
}

//return flux, given a time
double GRBobsmanager::flux(double time) const
{
  if(DEBUG) std::cout<<"flux at "<<time<<std::endl;
  double flux;

  if(m_grbdeleted) 
    flux = 0.0;
  else if(time <= m_startTime) 
    flux = 0.0;
  else 
    {
      if(!m_grbGenerated) const_cast<GRBobsmanager*>(this)->GenerateGRB();
      if(time < m_endTime)
	{
	  flux = m_spectrum->flux(time-m_startTime,m_MinPhotonEnergy);
	}
      else 
	{
	  const_cast<GRBobsmanager*>(this)->DeleteGRB();
	  flux = 0.0;
	}
    }
  if(DEBUG) std::cout<<"flux("<<time<<") = "<<flux<<std::endl;
  return flux;
}

double GRBobsmanager::interval(double time)
{  
  if(DEBUG) std::cout<<"interval at "<<time<<std::endl;
  double inte;
  if(m_grbdeleted) inte = 1e10;
  else if(time < m_startTime) 
    inte = m_startTime - time;
  else 
    {
      if(!m_grbGenerated) GenerateGRB();
      if(time == m_startTime) 
	{
	  inte = m_startTime - time + m_spectrum->interval(0.0,m_MinPhotonEnergy);
	}
      else if(time<m_endTime)
	{
	  inte = m_spectrum->interval(time - m_startTime,m_MinPhotonEnergy);
	}
      else 
	{
	  inte = 1e10;
	  DeleteGRB();
	}
    }
  if(DEBUG) std::cout<<"interval("<<time<<") = "<<inte<<std::endl;
  return inte;
}

double GRBobsmanager::energy(double time)
{
  if(DEBUG) std::cout<<"energy at "<<time<<std::endl;
  double ene=0.1;
  if(m_grbGenerated && !m_grbdeleted)
    ene = m_spectrum->energy(time-m_startTime,m_MinPhotonEnergy)*1.0e-3; //MeV
  if(DEBUG) std::cout<<"energy("<<time<<") = "<<ene<<std::endl;
  return ene;
}

double GRBobsmanager::parseParamList(std::string input, unsigned int index)
{
  std::vector<double> output;
  unsigned int i=0;
  for(;!input.empty() && i!=std::string::npos;){
    double f = ::atof( input.c_str() );
    output.push_back(f);
    i=input.find_first_of(",");
    input= input.substr(i+1);
  } 
  if(index>=(int)output.size()) return 0.0;
  return output[index];
}


