#include <iterator>
#include <map>
#include <vector>
#include <math.h>
#ifndef SpectObj_H
#define SpectObj_H

class SpectObj
{
 public:
  SpectObj(){;}
  SpectObj(double,double,double);
  SpectObj(std::map<double,double>::const_iterator,
	   std::map<double,double>::const_iterator);
  ~SpectObj(){;}
  
  double integrated_E_Flux(double emin, double emax);
  double integrated_Flux(double emin, double emax);
  
  double binSize(std::map<double,double>::const_iterator);
  
  inline double size()  { return m_spectrum.size(); }
  inline void   clear() { m_spectrum.clear();       }
  inline double getBinContent(double index) 
    {
      std::map<double,double>::const_iterator it = m_spectrum.begin();
      std::advance(it,index);
      return it->second;
    }
  
  inline double getBinValue(double index) 
    {
      std::map<double,double>::const_iterator it = m_spectrum.begin();
      std::advance(it,index);
      return it->first;
    }
  void clearSpectrumVector();
  //////////////////////////////////////////////////
  SpectObj operator+(const SpectObj inputObj);
  inline SpectObj operator+=(const SpectObj inputObj)
  {
    return (*this)+inputObj;
  }
  SpectObj operator*(const double value);
  inline SpectObj operator*=(const double value)
  {
    return (*this)*value;
  }
  SpectObj operator/(const double value);
  inline SpectObj operator/=(const double value)
  {
    return (*this)/value;
  }
  //////////////////////////////////////////////////
  std::vector<double> getEnergyVector(const double value = 1.);
  std::vector<double> getBinVector(const double value = 1.);
  std::vector<double> getSpectrumVector();
  
  SpectObj extractSub(std::map<double,double>::const_iterator,
                      std::map<double,double>::const_iterator);
  SpectObj SetSpectrum(double energy, double spectrum);
  SpectObj SetSpectrum(std::vector<double> energy, std::vector<double> spectrum);

  SpectObj AddSpectrum(SpectObj);
 private:
  std::map<double,double> m_spectrum;
  
};
#endif
