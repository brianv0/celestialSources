#ifndef TF1SPECTRUM_H
#define TF1SPECTRUM_H

#include "TF1.h"
#include "flux/Spectrum.h"
#include<map>

class TF1Spectrum : public Spectrum
{
 public:
  TF1Spectrum(const std::string& /*params*/);
  ~TF1Spectrum() {;}//{delete p_tf1;}
  
  float operator()(float )const 
    {
      return p_tf1.GetRandom();
    }

  std::string title() const 
    {
      return "Class creating a spectrum based on a formula";
    }

  const char * particleName() const
    {
      return m_particle_name.c_str();
    }

  
  double energy(double time)
    {
      return (*this)(time);
    }

 private:
   mutable TF1  p_tf1;
  std::map<std::string,std::string> m_parmap;
};

#endif
