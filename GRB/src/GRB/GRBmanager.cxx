#include "GRBmanager.h"
#include <iostream>

GRBmanager::GRBmanager(const std::string& params)
: m_params(params)
{
  m_FirstTime   = parseParamList(params,0);
  m_timeToWait  = parseParamList(params,1);
  m_initialTime = m_FirstTime;
  m_GRB = new GRBSpectrum("temp.txt");
  m_endTime = m_FirstTime+m_GRB->m_grbsim->Tmax();
  std::cout<<m_FirstTime<<std::endl;
  std::cout<<m_timeToWait<<std::endl;
  std::cout<<m_initialTime<<std::endl;
  std::cout<<m_endTime<<std::endl;
}

GRBmanager::~GRBmanager() 
{  
  std::cout<<"**************************************************"<<std::endl;
  delete m_GRB;
}

double GRBmanager::solidAngle() const
{
  return 1.0;
}

///return flux, given a time
double GRBmanager::flux(double time) const
{
  double flux;	  
  if(time < m_initialTime || (time > m_endTime)) flux = 1.;
 
  else flux = m_GRB->flux(time - m_initialTime);  
  std::cout<<"-----------------------------------"<<std::endl; 
  std::cout<<" flux @ internal time "<<time - m_initialTime<<" = "<<flux<<std::endl;
  std::cout<<"	time = "<<time<<std::endl;
  std::cout<<"	Initial time = "<<m_initialTime<<std::endl;
  std::cout<<"-----------------------------------"<<std::endl;
  return flux;
}

double GRBmanager::interval(double time)// const
{  
 
  double inte;
  	
  if(time < m_initialTime) inte = m_initialTime - time;
  else if (time<m_endTime)
	{
	  inte = m_GRB->interval(time - m_initialTime);
	}
  else if (time < m_endTime + m_timeToWait)
	{
	  delete m_GRB;
	  std::cout<<" NEW GRB "<<std::endl;
	  inte = m_timeToWait -(time - m_endTime);
	  m_GRB = new GRBSpectrum("temp.txt");
	  m_initialTime=time+inte;
	  m_endTime=m_initialTime+m_GRB->m_grbsim->Tmax();
//	  inte =  m_GRB->interval(0.0);
	}
  else
  	{
	 delete m_GRB;
	 std::cout<<" WARNING UNEXPECTED NEW GRB "<<std::endl;
	 inte = 0.0;
	 m_GRB = new GRBSpectrum("temp.txt");
	 m_initialTime=time+inte;
	 m_endTime=m_initialTime+m_GRB->m_grbsim->Tmax();
	} 	
  std::cout<<"-----------------------------------"<<std::endl;
  std::cout<<" Interval @ int time "<<time - m_initialTime<<std::endl;
  std::cout<<" 		time = "<<time<<std::endl;
  std::cout<<"		Initial Time = "<<m_initialTime<<std::endl;
  std::cout<<"		End Time = "<<m_endTime<<std::endl;
  std::cout<<" INTERVAL = "<<inte<<std::endl;
  std::cout<<"-----------------------------------"<<std::endl;
  return inte;
}

double GRBmanager::energySrc(HepRandomEngine* engine, double time)
{
  if(time < m_initialTime) return 0;
  return m_GRB->energySrc(engine,time - m_initialTime);
}

float GRBmanager::operator() (float u) const
{
  return (*m_GRB)(u);
}


float GRBmanager::parseParamList(std::string input, int index)
{
  std::vector<float> output;
  unsigned int i=0;
  for(;!input.empty() && i!=std::string::npos;){
    float f = ::atof( input.c_str() );
    output.push_back(f);
    i=input.find_first_of(",");
    input= input.substr(i+1);
  } 
  return output[index];
}


