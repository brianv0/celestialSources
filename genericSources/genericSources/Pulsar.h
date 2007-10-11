/**
 * @file Pulsar.h
 * @brief A pulsar with period and period derivative and light curve given
 * by a template file.
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef fluxSources_Pulsar_h
#define fluxSources_Pulsar_h

#include "flux/Spectrum.h"
#include "PeriodicSource.h"

/**
 * @class Pulsar
 * @brief Simulates a pulsar with period and period derivative and
 * light curve given by a template file.
 *
 * @author J. Chiang
 *
 * $Header$
 */

class Pulsar : public PeriodicSource {

public:

   Pulsar(const std::string &params);

   virtual ~Pulsar() {}

   /// @return Title describing the spectrum.
   virtual std::string title() const {return "Pulsar";}

   /// @return Interval to the next event (seconds)
   virtual double interval(double time);

   /// A useful static function for other classes (e.g., TransientTemplate)
   static void readLightCurve(const std::string & templateFile,
                              std::vector<double> & light_curve);

private:

   double m_meanFlux;
   double m_period;
   double m_pdot;
   double m_t0;
   double m_phi0;

   std::vector<double> m_phases;
   std::vector<double> m_integralDist;

   void computeIntegralDist(std::string templateFile);
   double period(double time) const;
   double residualTime(double time) const;
   double drawTime() const;
};

#endif // fluxSources_Pulsar_h
