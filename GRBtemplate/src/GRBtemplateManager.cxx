#include "GRBtemplate/GRBtemplateManager.h"
#include <iostream>
#include <fstream>
#include "flux/SpectrumFactory.h" 
#include "astro/EarthCoordinate.h"
#include "astro/GPS.h"
#include "astro/SkyDir.h"
#include "astro/PointingTransform.h"
#include "astro/JulianDate.h"

#define DEBUG 0 

ISpectrumFactory &GRBtemplateManagerFactory() 
{
  static SpectrumFactory<GRBtemplateManager> myFactory;
  return myFactory;
}


GRBtemplateManager::GRBtemplateManager(const std::string& params)
  : m_params(params)
{
  using astro::GPS;

  m_rnd = new TRandom();
  m_GenerateGBMOutputs = false;
  m_InputFileName = parseParamList(params,0);
  
  facilities::Util::expandEnvVar(&m_InputFileName);  
  
  m_startTime       = ::atof(parseParamList(params,1).c_str())+Spectrum::startTime();
  m_MinPhotonEnergy = ::atof(parseParamList(params,2).c_str())*1.0e3; //MeV

  if(::atof(parseParamList(params,3).c_str())!=0) m_GenerateGBMOutputs = true;
  m_theta_fow           = ::atof(parseParamList(params,4).c_str());
  
  
  
  //  m_par = new GRBtemplateParameters();
  m_GRBnumber = (long) floor(65540+m_startTime);
  //  facilities::Util::expandEnvVar(&m_InputFileName);
  
  m_theta = -1000000.0;
  
  //if(FOV>-90) std::cout<<"WARNING!! GRBtemplate: FOV = "<<FOV<<std::endl;
  m_theta_fow = 90.0 - m_theta_fow;
  while(m_theta < m_theta_fow - 2.0 || m_theta > m_theta_fow + 2.0)
    {
      double r1 = m_rnd->Uniform();
      double r2 = m_rnd->Uniform();
      m_l = 180.-360.*r1;
      m_b = ((180.0/M_PI)*acos(1.0-2.0*r2)-90.0);
      astro::SkyDir sky(m_l,m_b,astro::SkyDir::GALACTIC);
      HepVector3D skydir=sky.dir();
      try{
	HepRotation rottoglast = GPS::instance()->transformToGlast(m_startTime,GPS::CELESTIAL);
	HepVector3D scdir = rottoglast * skydir;
	m_ra    = sky.ra();
	m_dec   = sky.dec();
	double zenithCosTheta=cos(scdir.theta());
	m_theta = 90. - scdir.theta()*180.0/M_PI; // theta=0 -> XY plane, theta=90 -> Z
	m_phi   = scdir.phi()*180.0/M_PI;
	m_grbdeleted      = false;
	m_grbocculted = (zenithCosTheta < -0.4); // this is hardcoded  in FluxSource.cxx
      }
      catch(const std::exception &e)
	{
	  m_grbdeleted=true;
	  break;
	}
      m_GalDir=std::make_pair(m_l,m_b);
    }
  m_grbGenerated    = false;
  //////////////////////////////////////////////////
}

std::string GRBtemplateManager::GetGRBname()
{
  ////DATE AND GRBNAME
  /// GRB are named as customary GRBOBSYYMMDDXXXX
  //  The mission start and launch date are retrieved, 
  //  and the current burst's start time is added to them to form a Julian Date.
  //  Note: the argument of the JulianDate constructor is in day.

  astro::JulianDate JD(astro::JulianDate::missionStart() +
		       m_startTime/astro::JulianDate::secondsPerDay);

  int An,Me,Gio;
  m_Rest=((int) m_startTime % 86400) + m_startTime - floor(m_startTime);
  m_Frac=(m_Rest/86400.);
  int FracI=(int)(m_Frac*1000.0);
  double utc;
  JD.getGregorianDate(An,Me,Gio,utc);  
  An-=2000;

  std::ostringstream ostr;
  
  if(An<10) 
    {
      ostr<<"0";
      ostr<<An;
    }
  else
    {
      ostr<<An;
    }
  
  if(Me<10) 
    {
      ostr<<"0";
      ostr<<Me;
    }
  else
    {
      ostr<<Me;
    }
  
  if(Gio<10) 
    {
      ostr<<"0";
      ostr<<Gio;
    }
  else
    {
      ostr<<Gio;
    }
  
  if (FracI<10)
    {
      ostr<<"00";
      ostr<<FracI;
    }
  else if(FracI<100) 
    {
      ostr<<"0";
      ostr<<FracI;
    }
  else 
    ostr<<FracI;

  std::string GRBname = ostr.str();
  
  if(DEBUG) std::cout<<"GENERATE GRB TEMPLATE ("<<GRBname<<")"<<std::endl;
  
  //..................................................  //
  return GRBname;
}


GRBtemplateManager::~GRBtemplateManager() 
{  
  delete m_rnd;
  DeleteGRB();
}

void GRBtemplateManager::GenerateGRB()
{
  if(m_grbdeleted) return;
  /////GRB GRNERATION////////////////////////////////
  m_GRB      = new  GRBtemplateSim(m_InputFileName);
  m_spectrum = new  SpectObj(m_GRB->MakeGRB(),0);
  m_spectrum->SetAreaDetector(EventSource::totalArea());
  //////////////////////////////////////////////////
  m_endTime   = m_startTime + m_GRB->Tmax();
  std::string GRBname = GetGRBname();
  
  if(m_grbocculted)
    {
      std::cout<<"This GRB is occulted by the Earth"<<std::endl;
    }
  else if(m_GenerateGBMOutputs)
    {
      m_GRB->SaveGBMDefinition(GRBname,m_ra,m_dec,m_theta,m_phi,m_Rest);
      m_GRB->GetGBMFlux(GRBname);
    }
  
  TString name = "GRBTMP_";
  name+=GRBname; 
  name+="_PAR.txt";
  
  /*
    std::ofstream os(name,std::ios::out);
    os<<m_startTime<<" "<<m_endTime<<" "<<m_l<<" "<<m_b<<" "<<m_theta<<" "<<m_phi<<" "<<m_fluence<<" "<<" "<<m_alpha<<" "<<m_beta<<std::endl;
    os.close();
  */
  std::cout<<"Template GRB"<<GRBname<<" t start "<<m_startTime<<", tend "<<m_endTime
	   <<" l,b = "<<m_l<<", "<<m_b<<" elevation,phi(deg) = "<<m_theta<<", "<<m_phi<<std::endl;
  m_grbGenerated=true;
}

void GRBtemplateManager::DeleteGRB()
{
  if(m_grbGenerated)
    {
      delete m_GRB;
      delete m_spectrum;
    }
  m_grbdeleted=true;
}

//return flux, given a time
double GRBtemplateManager::flux(double time) const
{
  if(DEBUG) std::cout<<"Flux at "<<time<<" "<<m_startTime<<" "<<m_grbGenerated<<std::endl;
  double flux= 0.0;
  
  if(m_grbdeleted) 
    return flux;
  else if(time <= m_startTime || time >  m_startTime+5000) 
    {
      if(m_grbGenerated) const_cast<GRBtemplateManager*>(this)->DeleteGRB();
      return flux;
    }
  else 
    {
      if(!m_grbGenerated) const_cast<GRBtemplateManager*>(this)->GenerateGRB();
      flux = m_spectrum->flux(time-m_startTime,m_MinPhotonEnergy);
    }
  if(DEBUG) std::cout<<"flux("<<time<<") = "<<flux<<std::endl;
  return flux;
}

double GRBtemplateManager::interval(double time)
{  
  if(DEBUG) std::cout<<"interval at "<<time<<std::endl;
  double inte;
  if(m_grbdeleted) 
    {
      inte = 1e10;
    }
  else if(time < m_startTime)
    { 
      inte = m_startTime - time;
    }
  else if(time<m_startTime + 5000) //Confidential limit
    {
      if(!m_grbGenerated) const_cast<GRBtemplateManager*>(this)->GenerateGRB();
      inte = m_spectrum->interval(time - m_startTime,m_MinPhotonEnergy);
    }
  else  
    {
      inte = 1e10;
      if(m_grbGenerated) const_cast<GRBtemplateManager*>(this)->DeleteGRB();
    }
  
  if(DEBUG) std::cout<<"Interval("<<time<<") = "<<inte<<std::endl;
  return inte;
}

double GRBtemplateManager::energy(double time)
{
  if(DEBUG) std::cout<<"energy at "<<time<<std::endl;
  double ene=0.1;
  if(m_grbGenerated && !m_grbdeleted)
    ene = m_spectrum->energy(time-m_startTime,m_MinPhotonEnergy)*1.0e-3; //MeV
  if(DEBUG) std::cout<<"energy("<<time<<") = "<<ene<<std::endl;
  return ene;
}

std::string GRBtemplateManager::parseParamList(std::string input, unsigned int index)
{
  std::vector<std::string> output;
  unsigned int i=0;
  for(;!input.empty() && i!=std::string::npos;){
    
    i=input.find_first_of(",");
    std::string f = input.substr(0,i);
    output.push_back(f);
    input= input.substr(i+1);
  } 
  if(index>=(unsigned int)output.size()) return "0.0";
  return output[index];
}


