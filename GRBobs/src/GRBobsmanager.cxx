#include "GRBobs/GRBobsmanager.h"
#include <iostream>
#include <fstream>
#include "flux/SpectrumFactory.h" 

ISpectrumFactory &GRBobsmanagerFactory() 
{
  static SpectrumFactory<GRBobsmanager> myFactory;
  return myFactory;
}


GRBobsmanager::GRBobsmanager(const std::string& params)
  : m_params(params)
{
  
  facilities::Util::expandEnvVar(&paramFile);
  

  m_l               = parseParamList(params,0);
  m_b               = parseParamList(params,1);

  m_startTime       = parseParamList(params,2);
  m_fluence         = parseParamList(params,3);
  m_Npulses         = (int) parseParamList(params,4);
  m_alpha           = parseParamList(params,5);
  m_beta            = parseParamList(params,6);
  m_MinPhotonEnergy = parseParamList(params,7)*1.0e3; //MeV
  
  m_par = new GRBobsParameters();
  m_GRBnumber = 65540+m_Npulses;
  m_par->SetGRBNumber(65540+m_Npulses);  
  m_par->SetNumberOfPulses(m_Npulses);
  m_par->SetAlphaBeta(m_alpha,m_beta);
  m_par->SetMinPhotonEnergy(m_MinPhotonEnergy); //keV
  m_par->SetFluence(m_fluence);

  m_par->SetGalDir(m_l,m_b);

  m_grbGenerated    = false;
  //////////////////////////////////////////////////
  /*
    m_par = new GRBobsParameters();
  m_GRBnumber = 65540+m_Npulses;
  m_par->SetGRBNumber(65540+m_Npulses);  
  m_par->SetNumberOfPulses(m_Npulses);
  m_par->SetAlphaBeta(alpha,beta);
  m_par->SetMinPhotonEnergy(m_MinPhotonEnergy); //keV
  m_par->SetFluence(m_fluence);
  
  //////////////////////////////////////////////////
  m_GRB      = new  GRBobsSim(m_par);
  m_spectrum = new  SpectObj(m_GRB->MakeGRB(),0);
  m_spectrum->SetAreaDetector(EventSource::totalArea());
  //////////////////////////////////////////////////
  m_endTime   = m_startTime + m_GRB->Tmax();

  std::ofstream os("grbobs_generated.txt",std::ios::app);
  os<<m_GRBnumber<<" "<<m_startTime<<" "<<m_endTime<<" "<<m_fluence<<" "<<m_Npulses<<" "<<alpha<<" "<<beta<<std::endl;
  std::cout<<" Generate GRB obs ("<<m_GRBnumber<<"), ts = "<<m_startTime<<" te ="<<m_endTime<<" F= "<<m_fluence<<" N = "<<m_Npulses<<" a= "<<alpha<<" b= "<<beta<<std::endl;
  */
}

GRBobsmanager::~GRBobsmanager() 
{  
  delete m_par;
  delete m_GRB;
  delete m_spectrum;
}

void GRBobsmanager::GenerateGRB()
{

  //////////////////////////////////////////////////
  m_GRB      = new  GRBobsSim(m_par);
  m_spectrum = new  SpectObj(m_GRB->MakeGRB(),0);
  m_spectrum->SetAreaDetector(EventSource::totalArea());
  //////////////////////////////////////////////////
  m_endTime   = m_startTime + m_GRB->Tmax();
  
  std::ofstream os("grbobs_generated.txt",std::ios::app);
  os<<m_GRBnumber<<" "<<m_startTime<<" "<<m_endTime<<" "<<m_fluence<<" "<<m_Npulses<<" "<<m_alpha<<" "<<m_beta<<std::endl;
  std::cout<<" Generate GRB obs ("<<m_GRBnumber<<"), l = "<<m_l<<" b = "<<m_b<<" ts = "<<m_startTime<<" te ="<<m_endTime<<" F= "<<m_fluence<<" N = "<<m_Npulses<<" a= "<<m_alpha<<" b= "<<m_beta<<std::endl;
  m_grbGenerated=true;
	    
}

//return flux, given a time
double GRBobsmanager::flux(double time) const
{
  double flux;
  if(time <= m_startTime) 
    {
      flux = 0.0;
    }
  else 
    {
      if(!m_grbGenerated) GenerateGRB();
      flux = m_spectrum->flux(time-m_startTime,m_MinPhotonEnergy);
    }
  //  std::cout<<"flux ("<<time<<")= "<<flux<<std::endl;
  return flux;
}

double GRBobsmanager::interval(double time)
{  

  double inte;
  if(time < m_startTime) 
    inte = m_startTime - time;
  else 
    {
      if(!m_grbGenerated) GenerateGRB();
      if(time <= m_startTime) 
	{
	  inte = m_startTime - time + m_spectrum->interval(0.0,m_MinPhotonEnergy);
	}
      else if(time<m_endTime)
	{
	  if(! m_grbGenerated) GenerateGRB();
	  inte = m_spectrum->interval(time - m_startTime,m_MinPhotonEnergy);
	}
      else 
	{
	  delete m_par;
	  delete m_GRB;
	  delete m_spectrum;
	  inte = 1e10;
	}
    }
  //  std::cout<<"inte= "<< inte<<std::endl;
  return inte;
}

double GRBobsmanager::energy(double time)
{
  double ene;
  if(!m_grbGenerated)
    ene = 0.0;//m_MinPhotonEnergy * 1.0e-3; //MeV;
  else 
    ene = m_spectrum->energy(time-m_startTime,m_MinPhotonEnergy)*1.0e-3; //MeV
  std::cout<<"ene("<<time<<")= "<<ene<<std::endl;
  return ene;
}

double GRBobsmanager::parseParamList(std::string input, int index)
{
  std::vector<double> output;
  unsigned int i=0;
  for(;!input.empty() && i!=std::string::npos;){
    double f = ::atof( input.c_str() );
    output.push_back(f);
    i=input.find_first_of(",");
    input= input.substr(i+1);
  } 
  if(index>=output.size()) return 0.0;
  return output[index];
}


