#include <iostream>
#include "Pulsar/PulsarSpectrum.h"
#include "SpectObj/SpectObj.h"
#include "flux/SpectrumFactory.h"

ISpectrumFactory &PulsarSpectrumFactory() 
 {
   static SpectrumFactory<PulsarSpectrum> myFactory;
   return myFactory;
 }

PulsarSpectrum::PulsarSpectrum(const std::string& params)
  : m_params(params)
{
  std::cout<<m_fluence<<std::endl;
  m_fluence = parseParamList(params,1);
  m_period  = parseParamList(params,2);

  double ppar1 = 1e6;
  double ppar2 = 8e6;
  double ppar3 = -1.62;
  double ppar4 = 1.7;
  
  m_Pulsar   = new PulsarSim(m_fluence,m_period,2);
  m_spectrum = new SpectObj(m_Pulsar->PSRPolarCapPhen(ppar1,ppar2,ppar3,ppar4),1);
  //////////////////////////////////////////////////
  
}

PulsarSpectrum::~PulsarSpectrum() 
{  
  delete m_Pulsar;
  delete m_spectrum;
}

//return flux, given a time
double PulsarSpectrum::flux(double time) const
{
  double InternalTime = TMath::ASin(sin(2*TMath::Pi()*time/m_period))/(2*TMath::Pi());
  double flux;	  
  flux = m_spectrum->flux(InternalTime,cst::enph);
  return flux;
}

double PulsarSpectrum::interval(double time)
{  
  double inte;  
  double InternalTime = TMath::ASin(sin(2*TMath::Pi()*time/m_period))/(2*TMath::Pi());
  
  inte = m_spectrum->interval(InternalTime,cst::enph);
  return inte;
}

double PulsarSpectrum::energy(double time)
{
  double InternalTime = TMath::ASin(sin(2*TMath::Pi()*time/m_period))/(2*TMath::Pi());
  return m_spectrum->energy(InternalTime,cst::enph)*1.0e-3; //MeV
}

double PulsarSpectrum::parseParamList(std::string input, int index)
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


