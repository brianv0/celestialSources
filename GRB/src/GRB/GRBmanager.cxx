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
  cout<<m_FirstTime<<endl;
  cout<<m_timeToWait<<endl;
  cout<<m_initialTime<<endl;
  cout<<m_endTime<<endl;
}

GRBmanager::~GRBmanager() 
{  
  cout<<"**************************************************"<<endl;
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
  cout<<"-----------------------------------"<<endl; 
  cout<<" flux @ internal time "<<time - m_initialTime<<" = "<<flux<<endl;
  cout<<"	time = "<<time<<endl;
  cout<<"	Initial time = "<<m_initialTime<<endl;
  cout<<"-----------------------------------"<<endl;
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
	  cout<<" NEW GRB "<<endl;
	  inte = m_timeToWait -(time - m_endTime);
	  m_GRB = new GRBSpectrum("temp.txt");
	  m_initialTime=time+inte;
	  m_endTime=m_initialTime+m_GRB->m_grbsim->Tmax();
//	  inte =  m_GRB->interval(0.0);
	}
  else
  	{
	 delete m_GRB;
	 cout<<" WARNING UNEXPECTED NEW GRB "<<endl;
	 inte = 0.0;
	 m_GRB = new GRBSpectrum("temp.txt");
	 m_initialTime=time+inte;
	 m_endTime=m_initialTime+m_GRB->m_grbsim->Tmax();
	} 	
  cout<<"-----------------------------------"<<endl;
  cout<<" Interval @ int time "<<time - m_initialTime<<endl;
  cout<<" 		time = "<<time<<endl;
  cout<<"		Initial Time = "<<m_initialTime<<endl;
  cout<<"		End Time = "<<m_endTime<<endl;
  cout<<" INTERVAL = "<<inte<<endl;
  cout<<"-----------------------------------"<<endl;
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


