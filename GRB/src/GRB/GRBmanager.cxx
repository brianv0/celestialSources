#include "GRBmanager.h"
#include <iostream>
#include <fstream>
#include "flux/SpectrumFactory.h" 


ISpectrumFactory &GRBmanagerFactory() 
 {
   static SpectrumFactory<GRBmanager> myFactory;
   return myFactory;
 }
 
 
GRBmanager::GRBmanager(const std::string& params)
  : m_params(params)
{
  m_Nbursts=1;
  paramFile = "$(GRBROOT)/src/test/GRBParam.txt";
  facilities::Util::expandEnvVar(&paramFile);
  
  m_startTime   = TMath::Max(0.,parseParamList(params,0));
  m_timeToWait  = TMath::Max(0.,parseParamList(params,1));
  m_enph        = TMath::Max(0.,parseParamList(params,2))*1000.0; //keV
  m_par = new Parameters();
  //////////////////////////////////////////////////
  m_par->ComputeParametersFromFile(paramFile,1); 
  
  m_GRB      = new  GRBSim(m_par);
  m_spectrum = new  SpectObj(m_GRB->Fireball(),0);
  m_spectrum->SetAreaDetector(EventSource::totalArea());
  //////////////////////////////////////////////////
  m_endTime   = m_startTime + m_GRB->Tmax();
  m_nextBurst = m_endTime   + m_timeToWait;
  m_fluence   = m_GRB->GetFluence();
  m_GRBnumber = m_GRB->GetGRBNumber();
  m_GRBdir    =  m_GRB->GRBdir();

    std::ofstream os("grb_generated.txt",std::ios::out);

  os<<m_GRBnumber<<" "<<m_startTime<<" "<<m_endTime<<" "<<m_fluence<<" "<<m_GRBdir.first<<" "<<m_GRBdir.second<<std::endl;
  std::cout<<"GRB ("<<m_GRBnumber<<") starting at "<<m_startTime<<", ending at time: "<<m_endTime
	   <<" l,b ="<<m_GRBdir.first<<", "<<m_GRBdir.second<<" Fluence = "<<m_fluence<<std::endl;
  if(m_enph<=0) 
    m_enph=cst::enph;


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
      m_startTime = m_nextBurst;
      //////////////////////////////////////////////////
      m_Nbursts++;
      m_par->ComputeParametersFromFile(paramFile,m_Nbursts);
      m_GRB      = new  GRBSim(m_par);
      m_spectrum = new  SpectObj(m_GRB->Fireball(),0);
      m_spectrum->SetAreaDetector(EventSource::totalArea());
      //////////////////////////////////////////////////
      m_endTime   = m_startTime + m_GRB->Tmax();
      m_nextBurst = m_endTime   + m_timeToWait;      
      m_fluence   = m_GRB->GetFluence();
      m_GRBnumber = m_GRB->GetGRBNumber();
      m_GRBdir    =  m_GRB->GRBdir();

      inte = m_startTime-time + m_spectrum->interval(0.0,m_enph);
      std::cout<<"GRB ("<<m_GRBnumber<<") starting at "<<m_startTime<<", ending at time: "<<m_endTime
	       <<" l,b ="<<m_GRBdir.first<<", "<<m_GRBdir.second<<" Fluence = "<<m_fluence<<std::endl;
      std::ofstream os("grb_generated.txt",std::ios::app);
      os<<m_GRBnumber<<" "<<m_startTime<<" "<<m_endTime<<" "<<m_fluence<<" "<<m_GRBdir.first<<" "<<m_GRBdir.second<<std::endl;
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


