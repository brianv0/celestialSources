/**
 * @file IsotropicFileSpectrum.h
 * @brief Isotropic diffuse emission using an ascii file to specify
 * the spectral shape.
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef mySpectrum_IsotropicFileSpectrum_h
#define mySpectrum_IsotropicFileSpectrum_h

class FileSpectrum;

#include "flux/Spectrum.h"

/**
 * @class IsotropicFileSpectrum
 *
 * @brief IsotropicFileSpectrum diffuse emission.
 *
 */

class IsotropicFileSpectrum : public Spectrum {

public:

   /// This constructor is required and used in FluxSource for
   /// "SpectrumClass" sources.
   IsotropicFileSpectrum(const std::string &params);

   virtual ~IsotropicFileSpectrum(){}

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
   virtual std::string title() const {return "IsotropicFileSpectrum";}

   /// @return Interval to the next event (seconds)
   virtual double interval(double time);

   /// @return Photon energy (MeV).
   virtual double energy(double time);

   /// @return Photon direction in (l, b).
   virtual std::pair<double, double> dir(double energy);

private:

   FileSpectrum * m_fileSpectrum;

   double m_ra;
   double m_dec;
   double m_cos_thetamax;

};

#endif // mySpectrum_IsotropicFileSpectrum_h
