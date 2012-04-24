/**
 * @file GaussianSpectrum.h
 * @brief Point source with Gaussian spectrum
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef fluxSources_GaussianSpectrum_h
#define fluxSources_GaussianSpectrum_h

#include "flux/Spectrum.h"

/**
 * @class GaussianSpectrum
 * @brief Point source with a Gaussian spectrum
 *
 */

class GaussianSpectrum : public Spectrum {

public:

   /// This constructor is required and used in FluxSource for
   /// "SpectrumClass" sources.
   GaussianSpectrum(const std::string & params);

   virtual ~GaussianSpectrum(){}

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
   virtual std::string title() const {return "GaussianSpectrum";}

   /// @return Interval to the next event (seconds)
   virtual double interval(double time);

   /// @return Photon energy (MeV).
   virtual double energy(double time);

   /// @return Photon direction in (l, b).
   virtual std::pair<double, double> dir(double energy);

private:

   double m_flux;
   double m_mean;
   double m_sigma;

   double m_l;
   double m_b;
};

#endif // fluxSources_GaussianSpectrum_h
