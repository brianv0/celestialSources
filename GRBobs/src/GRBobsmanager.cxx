#include "GRBobs/GRBobsmanager.h"
#include <iostream>

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
  
  m_startTime       = parseParamList(params,0);
  m_fluence         = parseParamList(params,1);
  m_Npulses         = (int) parseParamList(params,2);
  double alpha      = parseParamList(params,3);
  double beta       = parseParamList(params,4);
  m_MinPhotonEnergy = parseParamList(params,5)*1.0e3; //MeV
  
  m_par = new GRBobsParameters();

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
  std::cout<<" GRB starting at time: "<<m_startTime<<" GRB ending at time: "<<m_endTime<<std::endl;
  std::cout<<" EventSource::totalArea()= "<<EventSource::totalArea()<<std::endl;
}

GRBobsmanager::~GRBobsmanager() 
{  
  delete m_par;
  delete m_GRB;
  delete m_spectrum;
}

//return flux, given a time
double GRBobsmanager::flux(double time) const
{
  double flux;	  
  if(time <= m_startTime || (time > m_endTime)) flux = 0.0;
  else flux = m_spectrum->flux(time-m_startTime,m_MinPhotonEnergy);
  return flux;
}

double GRBobsmanager::interval(double time)
{  
  double inte;  
  if(time <= m_startTime) inte = m_startTime - time + m_spectrum->interval(0.0,m_MinPhotonEnergy);
  else if (time<m_endTime)
    {
      inte = m_spectrum->interval(time - m_startTime,m_MinPhotonEnergy);
    }
  else 
    {
      
      inte = 1e10;
    }
  return inte;
}

double GRBobsmanager::energy(double time)
{
  double ene = m_spectrum->energy(time-m_startTime,m_MinPhotonEnergy)*1.0e-3; //MeV
  //  std::cout<<time<<" "<<ene<<std::endl;
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


