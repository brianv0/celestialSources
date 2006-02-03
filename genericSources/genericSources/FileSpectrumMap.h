#ifndef FILESPECTRUMMAP_H
#define FILESPECTRUMMAP_H

#include "genericSources/FileSpectrum.h"
#include "genericSources/MapSource.h"
#include<map>

class FileSpectrumMap : public MapSource
{
 public:
  FileSpectrumMap(const std::string& /*params*/);
  ~FileSpectrumMap() {delete m_filespectrum;}
  

  /// Overload of operator method
  float operator()(float time)const 
    {
      return (*m_filespectrum)(time);
    }

  std::string title() const 
    {
      return "Class creating a spectrum based on a formula";
    }

  const char * particleName() const    { return m_particle_name.c_str(); }

  

  /// Overload of flux method to ensure proper call to m_flux
  /// @return Total flux (photons/m^2).
  /// @param time Simulation time in seconds.
  virtual double flux(double ) const  { return m_flux; } 


 private:
  FileSpectrum*  m_filespectrum;
  std::map<std::string,std::string> m_parmap;
};

#endif


