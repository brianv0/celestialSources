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

double SpectObj::integrated_E_Flux(double enmin, double enmax) const
{  
  // m_spectrum  is in ph/s/MeV/m^2
  // integrated_E_Flux  is in eV/s/m^2
  
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

double SpectObj::integrated_Flux(double enmin, double enmax) const
{
  // m_spectrum  is in ph/s/MeV/m^2
  // integrated_E_Flux  is in ph/s/m^2

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


double SpectObj::binSize(std::map<double,double>::const_iterator it) const
{
  // if(*it = it.last()) return 0;
  double left_end  = (*it++).first;
  double right_end = (*it--).first;
  return (right_end - left_end);
}

std::vector<double> SpectObj::getEnergyVector(const double value) const
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

std::vector<double> SpectObj::getSpectrumVector(const double value) const
{
  std::vector<double> SpectrumVector;
  std::map<double,double>::const_iterator it = m_spectrum.begin();
  while(it != m_spectrum.end()) 
    {
      SpectrumVector.push_back(((*it++).second)*value);
    }
  return SpectrumVector;
}

std::vector<double> SpectObj::getBinVector(const double value) const
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

SpectObj SpectObj::SetSpectrum(int index, double spectrum)
  //purpose and method 
  //inputs 
  //outputs
  //caveat: needs a safeguard for out of range energy value 
{
  std::map<double,double>::iterator it = m_spectrum.begin();
  std::advance(it,index);
  (it->second) = spectrum;
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

const SpectObj SpectObj::operator+(const SpectObj inputObj)
{
  std::map<double,double>::iterator it = m_spectrum.begin();
  int index = 0;
  while(it != m_spectrum.end()) 
    {
      //m_spectrum[(*it++).first]+= inputObj.getSpectrum(index);
      (*it++).second += inputObj.getSpectrum(index);
      /*
	if (index < m_spectrum.size()-1)
      	{
	m_spectrum[(*it).first] += inputObj.getBinContent(index);
	}
      */
      //it++;
      index++;
    }
  return SpectObj();
}

const SpectObj SpectObj::operator*(const double value)
{
  std::map<double,double>::iterator it = m_spectrum.begin();
  while(it != m_spectrum.end())  
    {
      (*it++).second *= value;
    }
  return SpectObj();
}

const SpectObj SpectObj::operator/(const double value)
{
  std::map<double,double>::iterator it = m_spectrum.begin();
  while(it != m_spectrum.end())  
    {
      (*it++).second /= value;
    }
  return SpectObj();
}

double SpectObj::DrawPhotonFromSpectrum(float u, 
					double enmin) const
{
  // m_spectrum  is in ph/s/MeV/m^2
  // integrated_E_Flux  is in eV/s/m^2

  std::map<double,double>::const_iterator it = m_spectrum.begin();
  // With std:iterator
  // std::vector<double>::iterator it = spctrmVec.begin();  ;
  double value=0.;
  std::vector<double> Integral;
  std::vector<double>::iterator it_vec;
  int minbin=0;
  int pos;
  for(pos = 0; pos < m_spectrum.size()-1;pos++)
    {
      std::pair<double,double> bin = *it;
      double energy = bin.first;  //  eV
      double content= bin.second; //  ph/s/MeV/m^2
      if(energy >= enmin)
	{
	  value+=content*binSize(it);
	  Integral.push_back(value);
	}
      else
	minbin++;	
      it++;
    }  
  if(Integral.back() <=0 ) return 0.0;
  
  for(it_vec=Integral.begin();it_vec!=Integral.end();++it_vec) 
    { 
      (*it_vec) /= Integral.back(); 	//Normalizing to one
    }
  int nabove = Integral.size();
  int nbelow = 0;
  int ibin   = 0;
  while(nabove-nbelow > 1) 
    {
      int middle = (nabove+nbelow)/2;
      if (u == Integral[middle-1]) 
        {
          ibin = middle-1;
          break;
        }
      if (u  < Integral[middle-1]) nabove = middle;
      else                         nbelow = middle;
      ibin = nbelow-1;
    }
  //STEP4: retrurns the centered value of the energy at bin position
  // determined at STEP 3 
  
  double ph = getEnergy(ibin+minbin)+
    (getEnergy(ibin+minbin+1)-getEnergy(ibin+minbin))*
    (Integral[ibin+1] - u)/(Integral[ibin+1] - Integral[ibin]);
  
  return ph*1.0e-9; //returns value in GeV
}

SpectObj SpectObj::em2obs(double gamma, double angle)
{
  // energy  -> gamma*(1+beta*cos(angle)) * energy`
  // flux [in ph/s/MeV] -> ( gamma*(1+beta*cos(angle)))^2 [ph/s/MeV]
  double EnergyTransformation = getEnergyTransformation(gamma,angle);
  //  m_spectrum
  return (*this)*pow(EnergyTransformation,2.);
}

SpectObj SpectObj::obs2em(double gamma, double angle)
{
  // energy  -> gamma*(1+beta*cos(angle)) * energy`
  // flux [in ph/s/MeV] -> ( gamma*(1+beta*cos(angle)))^2 [ph/s/MeV]
  double EnergyTransformation = getEnergyTransformation(gamma,angle);
  //  m_spectrum*pow(1./EnergyTransformation,2.);
  return (*this)*pow(1./EnergyTransformation,2.);
}

  
//////////////////////////////////////////////////
