/**
 * @file RadialSource.h
 * @brief A source with radial extent and azimuthal symmetry.
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef fluxSources_RadialSource_h
#define fluxSources_RadialSource_h

#include <string>
#include <vector>

#include "flux/Spectrum.h"

class FileSpectrum;

/**
 * @class RadialSource
 *
 * @brief A source with radial extent and azimuthal symmetry.
 *
 */

class RadialSource : public Spectrum {

public:

   /// This constructor is required and used in FluxSource for
   /// "SpectrumClass" sources.
   RadialSource(const std::string &params);

   virtual ~RadialSource() {
      delete m_spectrum;
   }

   /// @return Particle energy in MeV.
   /// @param xi Uniform random deviate on the unit interval.
   virtual float operator() (float xi) const;

   /// @return Particle type, "gamma".
   virtual const char * particleName() const {return "gamma";}

   /// @return Total flux (photons/m^2).
   /// @param time Simulation time in seconds.
   virtual double flux(double time) const;

   /// @return "Effective" solid angle (sr).
   virtual double solidAngle() const;

   /// @return Title describing the spectrum.
   virtual std::string title() const {return "RadialSource";}

   /// @return Interval to the next event (seconds)
   virtual double interval(double time);

   /// @return Photon energy (MeV).
   virtual double energy(double time);

   /// @return Photon direction in (l, b).
   virtual std::pair<double, double> dir(double energy);

private:

   FileSpectrum * m_spectrum;

   std::vector<double> m_thetas;
   std::vector<double> m_radialDist;

   void fillRadialDist(const std::string & infile);
   double drawOffsetAngle() const;

   CLHEP::HepRotation m_rot;
   CLHEP::Hep3Vector m_rotatedSrcVec;
};

#endif // fluxSources_RadialSource_h
