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
  //  Field Of View for generating bursts degrees above the XY plane.
  const double FOV=-20;

  using astro::GPS;
  m_GenerateGBMOutputs = false;
  facilities::Util::expandEnvVar(&paramFile);  
  
  m_startTime       = parseParamList(params,0);
  double duration   = parseParamList(params,1);
  m_fluence         = parseParamList(params,2);
  m_z               = parseParamList(params,3);
  m_alpha           = parseParamList(params,4);
  m_beta            = parseParamList(params,5);
  m_MinPhotonEnergy = parseParamList(params,6)*1.0e3; //MeV

  if(parseParamList(params,7)!=0) m_GenerateGBMOutputs = true;

  m_LATphotons  = parseParamList(params,8);
  m_EC_delay    = parseParamList(params,9);
  m_EC_duration = parseParamList(params,10);
  m_Energy_CO   = parseParamList(params,11);
  
  m_par = new GRBobsParameters();

  m_GRBnumber = (long) floor(65540+m_startTime);

  m_par->SetGRBNumber(m_GRBnumber);
  m_par->SetFluence(m_fluence);
  m_par->SetPeakFlux(m_fluence);
  m_par->SetDuration(duration);
  m_par->SetAlphaBeta(m_alpha,m_beta);
  m_par->SetMinPhotonEnergy(m_MinPhotonEnergy); //keV  
  m_par->SetRedshift(m_z);
  m_theta = -100.0;
  
  double launch = Spectrum::startTime();

  while(m_theta<FOV)
    {
      m_par->SetGalDir(-200,-200); //this generates random direction in the sky
      m_GalDir = m_par->GetGalDir();      
      m_l = m_GalDir.first;
      m_b = m_GalDir.second;
      
      astro::SkyDir sky(m_l,m_b,astro::SkyDir::GALACTIC);
      HepVector3D skydir=sky.dir();
      HepRotation rottoglast = GPS::instance()->transformToGlast(m_startTime+launch,GPS::CELESTIAL);
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
  /// GRB are named as customary GRBYYMMDDXXXX
  //  The mission start and launch date are retrieved, 
  //  and the current burst's start time is added to them to form a Julian Date.
  //  Note: the argument of the JulianDate constructor is in day.

  astro::JulianDate JD(astro::JulianDate::missionStart() +
		       (Spectrum::startTime()+m_startTime)/astro::JulianDate::secondsPerDay);

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
  m_GRB      = new GRBobsSim(m_par);
  TH2D *h;
  h = m_GRB->MakeGRB();
  m_spectrum = new SpectObj(m_GRB->CutOff(h,m_Energy_CO) , 0, m_z);
  m_spectrum->SetAreaDetector(EventSource::totalArea());
  m_endTime    = m_startTime  + m_GRB->Tmax();
  
  PromptEmission = new SpectralComponent(m_spectrum,m_startTime,m_endTime);
  //std::cout<<"Generate PROMPT emission ("<<m_startTime<<" "<<m_endTime<<")"<<std::endl;
  if(m_LATphotons>0)
    {
      m_startTime_EC = m_startTime    + m_EC_delay;
      m_endTime_EC   = m_startTime_EC + m_EC_duration;
      h = m_GRB->MakeGRB_ExtraComponent(m_EC_duration,m_LATphotons);
      m_spectrum1 = new SpectObj(m_GRB->CutOff(h,m_Energy_CO),0, m_z);
      m_spectrum1->SetAreaDetector(EventSource::totalArea());
      AfterGlowEmission = new SpectralComponent(m_spectrum1,m_startTime_EC,m_endTime_EC);
      //      std::cout<<"Generate AFTERGLOW Emission ("<<m_startTime_EC<<" "<<m_endTime_EC<<")"<<std::endl;
    }
  
  //////////////////////////////////////////////////
  TString GRBname = GetGRBname();

  if(m_GenerateGBMOutputs)
    {
      m_GRB->SaveGBMDefinition(GRBname,m_ra,m_dec,m_theta,m_phi,m_Rest);
      m_GRB->GetGBMFlux(GRBname);
    }
  
  TString name = "GRBOBS_";
  name+=GRBname; 
  name+="_PAR.txt";
  //..................................................//
  std::ofstream os(name,std::ios::out);  
  os<<m_startTime<<" "<<m_endTime<<" "<<m_l<<" "<<m_b<<" "<<m_theta<<" "<<m_phi<<" "<<m_fluence<<" "<<m_alpha<<" "<<m_beta<<std::endl;
  os.close();
  
  std::cout<<"GRB"<<GRBname;
  if( m_z >0)  std::cout<<" redshift = "<<m_z;
  std::cout<<" t start "<<m_startTime<<", tend "<<m_endTime
	   <<" l,b = "<<m_l<<", "<<m_b<<" elevation,phi(deg) = "<<m_theta<<", "<<m_phi;
  if(m_par->GetNormType()=='P')
    std::cout<<" Peak Flux = "<<m_fluence<<" 1/cm^2/s "<<std::endl;
  else 
    std::cout<<" Fluence = "<<m_fluence<<" erg/cm^2"<<std::endl;
  if(m_LATphotons>0) 
    {
      m_endTime = TMath::Max(m_endTime, m_startTime + m_EC_delay  + m_EC_duration);
      std::cout<<"Generate AFTERGLOW Emission ("<<m_startTime_EC<<" "<<m_endTime_EC<<")"<<std::endl;
    }
  m_grbGenerated=true;
}

void GRBobsmanager::DeleteGRB()
{
  delete m_par;
  if(m_grbGenerated)
    {
      delete m_GRB;
      delete m_spectrum;
      delete PromptEmission;
      if(m_LATphotons>0) 
	{
	  delete AfterGlowEmission;
	  delete m_spectrum1;
	}
    }
  m_grbdeleted=true;
}

//return flux, given a time
double GRBobsmanager::flux(double time) const
{
  double flux = 0.0;
  //////////////////////////////////////////////////
  
  if(m_grbdeleted) 
    flux = 0.0;
  else if(time <= m_startTime) 
    flux = 0.0;
  else 
    {
      if(!m_grbGenerated) const_cast<GRBobsmanager*>(this)->GenerateGRB();
      if(time < m_endTime)
	{
	  flux   = PromptEmission->flux(time,m_MinPhotonEnergy);
	  if(m_LATphotons>0) flux  += AfterGlowEmission->flux(time,m_MinPhotonEnergy);
	}
      else 
	{
	  const_cast<GRBobsmanager*>(this)->DeleteGRB();
	  flux = 0.0;
	}
    }
  if(DEBUG && flux) std::cout<<"flux("<<time<<") = "<<flux<<std::endl;
  return flux;
}

double GRBobsmanager::interval(double time)
{  
  if(DEBUG) std::cout<<"interval at "<<time<<std::endl;
  double inte, inte_prompt, inte_ag;
  inte_ag=1e20;
  if(m_grbdeleted) inte = 1e10;
  else if(time < m_startTime) 
    inte = m_startTime - time;
  else 
    {
      if(!m_grbGenerated) GenerateGRB();
      if(time<m_endTime)
	{
	  inte_prompt = PromptEmission->interval(time,m_MinPhotonEnergy);
	  if(m_LATphotons>0) inte_ag = AfterGlowEmission->interval(time,m_MinPhotonEnergy);
	  if (inte_prompt<=inte_ag) 
	    { 
	      inte  = inte_prompt;
	      PromptEmission->SetResiduals(0.0);
	      if(m_LATphotons>0) AfterGlowEmission->SetResiduals(inte);
	    }
	  else
	    { 
	      inte  = inte_ag;
	      PromptEmission->SetResiduals(inte);
	      if(m_LATphotons>0) AfterGlowEmission->SetResiduals(0.0);
	    }
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
  double ene_prompt,ene_ag;
  ene_ag=0.0;
  if(m_grbGenerated && !m_grbdeleted)
    {
      if(m_LATphotons>0) ene_ag = AfterGlowEmission->energyMeV(time,m_MinPhotonEnergy);
      ene = TMath::Max(PromptEmission->energyMeV(time,m_MinPhotonEnergy),ene_ag);
    }
  if(DEBUG) 
    std::cout<<"energy("<<time<<") = "<<ene<<std::endl;
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


