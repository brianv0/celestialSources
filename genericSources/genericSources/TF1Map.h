#ifndef TF1MAP_H
#define TF1MAP_H

#include "TF1.h"
#include "genericSources/MapSource.h"
#include<map>

class TF1Map : public MapSource
{
 public:
  TF1Map(const std::string& /*params*/);
  ~TF1Map() {;}//{delete p_tf1;}
  

  /// Overload of operator method
  float operator()(float )const 
    {
      return p_tf1.GetRandom();
    }

  std::string title() const 
    {
      return "Class creating a spectrum based on a formula";
    }

  const char * particleName() const    {      return m_particle_name.c_str();    }

  
  double energy(double time)    {      return (*this)(time);    }


  /// Overload of flux method to ensure proper call to m_flux
  /// @return Total flux (photons/m^2).
  /// @param time Simulation time in seconds.
  virtual double flux(double ) const  { return m_flux; } 


 private:
  TF1  p_tf1;
  std::map<std::string,std::string> m_parmap;
};

#endif
