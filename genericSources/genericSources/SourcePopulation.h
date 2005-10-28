/**
 * @file SourcePopulation.h
 * @brief Draw from a population of steady point sources.
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef genericSources_SourcePopulation_h
#define genericSources_SourcePopulation_h

#include <string>
#include <utility>
#include <vector>

#include "flux/Spectrum.h"

namespace astro {
   class SkyDir;
}

namespace IRB {
   class EblAtten;
}

/**
 * @class SourcePopulation
 * @brief Draw from a population of steady point sources.  The sources
 * are read in from a standard flux-style xml file and so may be
 * modeled using a power-low or broken power-law.
 *
 * @author J. Chiang
 *
 * $Header$
 */

class SourcePopulation : public Spectrum {

public:

   SourcePopulation(const std::string & params);
   
   /// @return Particle energy in MeV.
   /// @param xi Uniform random deviate on the unit interval.
   virtual float operator() (float xi) const {
      (void)(xi);
      return m_currentEnergy;
   }

   /// @return Particle type, "gamma".
   virtual const char * particleName() const {
      return "gamma";
   }

   /// @return Title describing the spectrum.
   virtual std::string title() const {
      return "SourcePopulation";
   }

   /// @return Time interval to the next event.  This implementation
   /// returns -1 to instruct the callling class to use the cumulative
   /// flux value.  However, calling it causes the photon direction
   /// and energy to be set.
   virtual double interval(double time);

   /// @return Photon energy (MeV).
   virtual double energy(double time) {
      (void)(time);
      return m_currentEnergy;
   }

   /// @return photon direction in (l, b)
   /// @param energy This is the energy returned from a previous call
   /// to the energy(time) method (see the flux/ISpectrum
   /// declaration).  For this source, which has effectively a
   /// spatially varying spectra, this logic is *backwards*.  One
   /// would want to find the incident direction first, then find the
   /// energy.
   virtual std::pair<double, double> dir(double energy) {
      (void)(energy);
      return std::make_pair(m_l, m_b);
   }

private:

   double m_l;
   double m_b;
   double m_currentEnergy;

   std::vector<double> m_cumulativeFlux;

   class PointSource {
   public:
      PointSource(const astro::SkyDir & dir, double flux,
                  double gamma, double gamma2, double ebreak,
                  double emin, double emax, const IRB::EblAtten * tau=0);

      ~PointSource();

      /// @return Photon direction as a SkyDir object
      const astro::SkyDir & dir() const {
         return m_dir;
      }

      /// @return Draw the photon energy from the spectral model (MeV)
      double energy() const;

      /// @return Photon flux (#/m^2/s)
      double flux() const {
         return m_flux;
      }

   private:

      astro::SkyDir m_dir;
      double m_flux;
      double m_gamma;
      double m_gamma2;
      double m_ebreak;
      double m_emin;
      double m_emax;
      const IRB::EblAtten * m_tau;

      double m_part1;
      double m_part2;
      double m_frac;
   };

   std::vector<PointSource> m_sources;

   void readSourceFile(const std::string & input_file);

};

#endif // genericSources_SourcePopulation_h
