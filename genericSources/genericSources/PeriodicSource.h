/**
 * @file PeriodicSource.h
 * @brief Simple periodic source with sinusoidal modulation.
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef mySpectrum_PeriodicSource_h
#define mySpectrum_PeriodicSource_h

#include "flux/Spectrum.h"

//namespace fluxSources {

/**
 * @class PeriodicSource
 * @brief A periodic source with sinusoidal modulation for the flux
 * package.
 *
 * @author J. Chiang
 *
 * $Header$
 */

class PeriodicSource : public Spectrum {

public:

   /// This constructor is required and used in FluxSource for
   /// "SpectrumClass" sources.
   PeriodicSource(const std::string &params);

   virtual ~PeriodicSource(){}

   /// @return Particle energy in MeV.
   /// @param xi Uniform random deviate on the unit interval.
   virtual float operator() (float xi) const;

   /// @return Particle type, "gamma".
   virtual const char * particleName() const {return "gamma";}

   /// @return Title describing the spectrum.
   virtual std::string title() const {return "PeriodicSource";}

   /// @return Interval to the next event (seconds)
   virtual double interval(double time);

   /// @return Photon energy (MeV).
   virtual double energy(double time);

protected:

   PeriodicSource(double gamma, double emin, double emax) 
      : m_gamma(gamma), m_emin(emin), m_emax(emax) {}

   double m_gamma;
   double m_emin;
   double m_emax;

   void makeGrid(unsigned int n, double xmin, double xmax, 
                 std::vector<double> &x, bool makeLog=false);
   double interpolate(const std::vector<double> &x, 
                      const std::vector<double> &y,
                      double xx);

   /// Disable these virtual functions since they are not used by
   /// this source.
   virtual double flux(double) const {return 0;}
   virtual double solidAngle() const {return 0;}

   /// galactic_dir or celestial_dir tags will give the
   /// photon directions.
   virtual std::pair<double, double> dir(double) {
      return std::make_pair(0, 0);
   }

private:

   double m_flux;
   double m_period;
   double m_amplitude;
   double m_phi0;

   std::vector<double> m_arrTimes;
   std::vector<double> m_integralDist;

   void computeIntegralDistribution();

};

//} // namespace fluxSources

#endif // mySpectrum_PeriodicSource_h
