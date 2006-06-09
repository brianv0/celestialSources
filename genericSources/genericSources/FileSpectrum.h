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
   virtual float operator()(float x);
    
   virtual std::string title() const;

   virtual const char * particleName() const;

   inline const char * nameOf() const {
      return "FileSpectrum";
   }

private:

   typedef std::pair<double, double> EFpair_t;

   std::vector<EFpair_t> m_integralSpectrum;

   double read_file(const std::string & infile);

   double compute_integral_dist(const std::vector<double> & energies,
                                const std::vector<double> & fluxes);
   
};

#endif // genericeSources_FileSpectrum_H
