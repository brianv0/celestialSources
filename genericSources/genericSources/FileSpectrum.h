/** 
 * @file FileSpectrum.h

 * @brief FileSpectrum allows definition of a source based on a
 * spectrum saved in a file. The file is currently assumed to have two
 * columns: energy value and flux (normalized or not). Flux value can
 * be renormalized using the xml source description.  Energy in the
 * file can only be in MeV consistent with GLAST/LAT standards.
 *
 *  $Header$
 */

#ifndef genericSources_FileSpectrum_H
#define genericSources_FileSpectrum_H

#include <map>
#include <string>
#include <vector>

#include "flux/Spectrum.h"

/** 
 * \class FileSpectrum
 *
 * \brief Spectrum that reads its differential spectrum from a table
 * \author Theodore Hierath
 * 
 * $Header$
 */

class FileSpectrum : public Spectrum {

public:

   FileSpectrum(const std::string& params);
   
   /// @return total flux in #/m^2/s
   virtual double flux() const;

   /// @return total flux in #/m^2/s
   virtual double flux (double time) const;
    
   /// sample a single particle energy from the spectrum
   virtual float operator() (float);
    
   virtual std::string title() const;

   virtual const char * particleName() const;

   inline  const char * nameOf() const {return "FileSpectrum";}
   
private:

// normalization from the file: can be overwritten in the xml
   float m_fileflux;
   
   typedef std::pair<double, double> efpair_t;
   std::vector<efpair_t> m_integ_flux;

   std::vector<efpair_t> readFile(std::ifstream& input_file);

};

#endif // genericeSources_FileSpectrum_H
