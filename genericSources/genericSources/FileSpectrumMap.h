/**
 * @file FileSpectrumMap.h
 * @brief A source class that uses a FileSpectrum instance to define a spectrum via an ascii file,
 * and its localization in space via a Fits map, encapsulated in the MapSource base class
 * The class is instantiated based on the information parsed from
 * the params string, passed to the object vi the constructor. Here is an example of
 * such params string :
 * params="flux=17.,fitsFile=$(FLUXROOT)/sources/gas_gal.fits,specFile=$(GENERICSOURCESROOT)/data/dm120gev.dat,emin=100.,emax=1100,lonMin=-180,lonMax=180,latMin=-90,latMax=90"/>
 * emin,emax,lonMin,lonMax,latMin,latMax are used by MapSource.
 * @author Johann Cohen-Tanugi
 *
 * $Header$
 */

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
      return "FileSpectrumMap";
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


