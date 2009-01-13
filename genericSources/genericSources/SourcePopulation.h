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
 * modeled using a power-law or broken power-law.
 *
 * @author J. Chiang
 *
 * $Header$
 */

class SourcePopulation : public Spectrum {

public:

   SourcePopulation(const std::string & params);

   virtual ~SourcePopulation();
   
   /// @return Time interval to the next event
   virtual double interval(double time);

   /// @return Particle energy in MeV.
   /// @param xi Uniform random deviate on the unit interval.
   virtual float operator() (float xi);

   /// @return Particle type, "gamma".
   virtual const char * particleName() const {
      return "gamma";
   }

   /// @return Title describing the spectrum.
   virtual std::string title() const {
      return "SourcePopulation";
   }

   /// @return Photon energy (MeV).
   virtual double energy(double time);

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

   virtual std::string name() const;

private:

   IRB::EblAtten * m_tau;

   int m_idOffset;
   double m_l;
   double m_b;
   double m_currentEnergy;
   std::string m_name; /// name of current source

   std::vector<double> m_cumulativeFlux;

   class PointSource {
   public:
      PointSource(const astro::SkyDir & dir, double flux,
                  double gamma, double gamma2, double ebreak,
                  double emin, double emax, double zz);

      PointSource(const std::string & line);

      ~PointSource() {}

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

      static void setEblAtten(IRB::EblAtten * tau) {
         s_tau = tau;
      }

      double integral(double emin, double emax);

      std::string name()const{
          return m_name;
      }

   private:

      astro::SkyDir m_dir;
      double m_flux;
      double m_gamma;
      double m_gamma2;
      double m_ebreak;
      double m_emin;
      double m_emax;
      double m_z;

      static IRB::EblAtten * s_tau;

      double m_part1;
      double m_part2;
      double m_frac;
      std::string m_name; /// name of the source

      void setPowerLaw();

      double dnde(double energy) const;

      double attenuation(double energy) const;

      class DndeIntegrand {
      public:
         DndeIntegrand(const PointSource & pointSource) 
            : m_pointSource(pointSource) {}
         double operator()(double energy) const {
            return m_pointSource.dnde(energy);
         }
      private:
         const PointSource & m_pointSource;
      };

   };

   std::vector<PointSource> m_sources;

   void setEblAtten(const std::string & ebl_par);

   void readSourceFile(std::string input_file);

};

#endif // genericSources_SourcePopulation_h
