#include "SpectObj.h"

SpectObj::SpectObj(double const enmin, double const enmax, double const enstep) 
{
  double temp1=(enmax/enmin);
  double temp2=(1.0/(enstep));
  double denergy = pow(temp1,temp2);
  int en;
  for(en=0;en<=enstep;en++)
    {
      m_spectrum[enmin*pow(denergy,en)]=0.;
    }
}

SpectObj::SpectObj(std::map<double,double>::const_iterator start,
		   std::map<double,double>::const_iterator end)
{
  for(std::map<double,double>::const_iterator it = start; it!=end;it++)
    {
      m_spectrum[(*it).first] = (*it).second;
    }
}

double SpectObj::integrated_E_Flux(double enmin, double enmax)
{
  double value = 0.;
  std::map<double,double>::const_iterator it = m_spectrum.begin();
  while(it != m_spectrum.end())
    {
      
      std::pair<double,double> bin = *it;
      double energy = bin.first;
      double content= bin.second;
      if(energy>=enmin && energy<enmax)
        {
          value += energy*content*binSize(it);
        }
      it++;
    }
  value *= 1.0e-6;
  return value;
}


double SpectObj::integrated_Flux(double enmin, double enmax)
{
  
  double value = 0;
  std::map<double,double>::const_iterator it = m_spectrum.begin();
  while(it != m_spectrum.end())
    {
      std::pair<double,double> bin = *it;
      double energy = bin.first;
      double content= bin.second;
      if(energy>=enmin && energy<enmax)
        {
	  std::map<double,double>::const_iterator it2 = it;
	  it2++;
	  if(it2 == m_spectrum.end())
	    {
	      value = 0;
	    } 
	  else
	    {
	      value += content*binSize(it);
	    }
        }
      it++;
    }
  value *= 1.0e-6;
  return value;
}


double SpectObj::binSize(std::map<double,double>::const_iterator it)
{
  // if(*it = it.last()) return 0;
  double left_end  = (*it++).first;
  double right_end = (*it--).first;
  return (right_end - left_end);
}

std::vector<double> SpectObj::getEnergyVector(const double value)
{
  std::vector<double> EnergyVector;
  std::map<double,double>::const_iterator it = m_spectrum.begin();
  while(it != m_spectrum.end()) 
    {
      EnergyVector.push_back(((*it).first)*value);
      it++;
    }
  return EnergyVector;
}

std::vector<double> SpectObj::getSpectrumVector(const double value)
{
  std::vector<double> SpectrumVector;
  std::map<double,double>::const_iterator it = m_spectrum.begin();
  while(it != m_spectrum.end()) 
    {
      SpectrumVector.push_back(((*it++).second)*value);
    }
  return SpectrumVector;
}

std::vector<double> SpectObj::getBinVector(const double value)
{
  std::vector<double> BinVector;
  std::vector<double> EnergyVector = getEnergyVector(value);
  std::vector<double>::iterator it = EnergyVector.begin();
  while(it+1!= EnergyVector.end()) 
    {
      BinVector.push_back(*(it+1)-(*it));
      it++;
    }
  return BinVector;
}

SpectObj SpectObj::extractSub(std::map<double,double>::const_iterator start,
                              std::map<double,double>::const_iterator end)
{
  return SpectObj(start,end);
}

SpectObj SpectObj::SetSpectrum(double energy, double spectrum)
  //purpose and method 
  //inputs 
  //outputs
  //caveat: needs a safeguard for out of range energy value 
{
  m_spectrum[energy] = spectrum;
  return (*this);
}

SpectObj SpectObj::SetSpectrum(std::vector<double> energy, std::vector<double> spectrum)
  //purpose and method 
  //inputs 
  //outputs
  //caveat: needs a safeguard for out of range energy and spectrum values
{ 
  std::vector<double>::const_iterator ite = energy.begin();
  std::vector<double>::const_iterator its = spectrum.begin();
  while (ite != energy.end())
    {
      m_spectrum[(*ite++)] = (*its++);
    }
  return (*this);
}

//////////////////////////////////////////////////

SpectObj SpectObj::operator+(SpectObj inputObj)
{
  std::map<double,double>::const_iterator it = m_spectrum.begin();
  unsigned int index = 0;
  while(it != m_spectrum.end()) 
    {
      if (index < m_spectrum.size()-1)
	{
	  m_spectrum[(*it).first] += inputObj.getBinContent(index);
	}
      it++;
      index++;
    }
  return SpectObj();
}

SpectObj SpectObj::operator*(const double value)
{
  std::map<double,double>::iterator it = m_spectrum.begin();
  while(it != m_spectrum.end())  
    {
      (*it++).second *= value;
    }
  return SpectObj();
}

SpectObj SpectObj::operator/(const double value)
{
  std::map<double,double>::iterator it = m_spectrum.begin();
  while(it != m_spectrum.end())  
    {
      (*it++).second /= value;
    }
  return SpectObj();
}
//////////////////////////////////////////////////
