/**
 * @file SpectralTransient.h
 * @brief A flaring source whose light curve and spectral variability
 * is given by a template file.
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef genericSources_SpectralTransient_h
#define genericSources_SpectralTransient_h

#include "flux/Spectrum.h"

namespace IRB {
   class EblAtten;
}

/**
 * @class SpectralTransient
 *
 * @brief A flaring source whose light curve shape is given by a
 * template file.  The duration and mean photon flux are given as
 * parameters.  The spectrum during each interval defined in the
 * template file is given as a broken power-law.
 *
 * @author J. Chiang
 *
 * $Header$
 */

class SpectralTransient : public Spectrum {

public:

   SpectralTransient(const std::string & params);
   
   virtual ~SpectralTransient() {}

   /// @return Particle energy in MeV.
   /// @param xi Uniform random deviate on the unit interval.
   virtual float operator() (float xi) const {
      (void)(xi);
      return m_currentEnergy;
   }

   /// @return Particle type, "gamma".
   virtual const char * particleName() const {return "gamma";}

   /// @return Title describing the spectrum.
   virtual std::string title() const {return "SpectralTransient";}

   /// @return Interval to the next event (seconds)
   virtual double interval(double time);

   /// @return Photon energy (MeV).
   virtual double energy(double time) {
      (void)(time);
      return m_currentEnergy;
   }

   const std::vector<std::pair<double, double> > & events() const {
      return m_events;
   }

private:

   double m_flux;
   double m_tstart;
   double m_tstop;
   double m_emin;
   double m_emax;
   int m_lc;
   float m_z;

   IRB::EblAtten * m_tau;

   double m_currentEnergy;

   std::vector<std::pair<double, double> > m_events;

   void createEvents(std::string templateFile);

   class ModelInterval {
   public:

      ModelInterval() {}

      /// Read data members from a line in the template file.
      ModelInterval(const std::string & line, double emin, double emax);

      /// Pass the data members via an ordered vector.
      ModelInterval(const std::vector<double> & data,
                    double emin, double emax);

      /// Fractional start time of the interval; the entire light curve 
      /// will be rescaled to fit the interval [m_tstart, m_tstop]
      double startTime;
      /// Fraction stop time.
      double stopTime;
      /// Energy integrated photon flux (#/cm^2-s)
      double flux;
      /// Lower energy photon index
      double gamma1;
      /// Higher energy photon index
      double gamma2;
      /// Break energy (MeV)
      double ebreak;

      bool drawBelowBreak(double xi) const {
         return xi < m_lowerFraction;
      }

   private:
      double m_lowerFraction;

      void brokenPowerLawFractions(double emin, double emax);
   };

   std::vector<ModelInterval> m_lightCurve;

   void readLightCurve(const std::string & templateFile);
   void readFitsLightCurve(const std::string & templateFile);

   void rescaleLightCurve();

   std::pair<double, double> 
   drawEvent(const std::vector<double> & integralDist) const;

   double drawEnergy(const ModelInterval & interval) const;

   static bool compareEventTime(const std::pair<double, double> & x,
                                const std::pair<double, double> & y) {
      return x.first < y.first;
   }
};

#endif // genericSources_SpectralTransient_h
