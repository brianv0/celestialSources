/**
 * @file GaussianSource.h
 * @brief Test spectrum for use with flux package.
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef fluxSources_GaussianSource_h
#define fluxSources_GaussianSource_h

#include "flux/Spectrum.h"

/**
 * @class GaussianSource
 * @brief Simple class to exercise flux package.
 *
 * @author J. Chiang
 *
 * $Header$
 */

class GaussianSource : public Spectrum {

public:

   /// This constructor is required and used in FluxSource for
   /// "SpectrumClass" sources.
   GaussianSource(const std::string &params);

   virtual ~GaussianSource(){}

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
   virtual std::string title() const {return "GaussianSource";}

   /// @return Interval to the next event (seconds)
   virtual double interval(double time);

   /// @return Photon energy (MeV).
   virtual double energy(double time);

   /// @return Photon direction in (l, b).
   virtual std::pair<double, double> dir(double energy);

private:

   double m_flux;
   double m_gamma;
   double m_major;
   double m_minor;
   double m_posAngle;
   double m_emin;
   double m_emax;

   HepRotation m_rot;
   Hep3Vector m_rotatedSrcVec;

};

#endif // fluxSources_GaussianSource_h
