/** 
 * @file FileSpectrum.h
 * @brief FileSpectrum allows definition of a source based on a
 * spectrum from an ascii file. The file must have two columns: energy
 * and flux (normalized or not). The flux values can be renormalized
 * using the xml source description.  Energy values in the file must
 * be in MeV, consistent with GLAST/LAT standards.
 *
 *  $Header$
 */

#ifndef genericSources_FileSpectrum_H
#define genericSources_FileSpectrum_H

#include <deque>
#include <string>
#include <vector>

#include "flux/Spectrum.h"

/** 
 * @class FileSpectrum
 *
 * @brief Spectrum that reads its differential spectrum from a
 * table. This version has been massively refactored from the original
 * implementation in flux/FILESpectrum by Theodore Hierath.
 *
 */

class FileSpectrum : public Spectrum {

public:

   FileSpectrum(const std::string & params);
   
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

   double m_emin;
   double m_emax;
   std::deque<double> m_energies;
   std::deque<double> m_dnde;
   std::deque<double> m_integralSpectrum;

   double read_file(const std::string & infile);
   void reset_ebounds();
   double compute_integral_dist();
   
};

#endif // genericeSources_FileSpectrum_H
