/**
 * @file FitsTransient.h
 * @brief A flaring source whose spectral evolution is given by a FITS
 * binary table.
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef genericSources_FitsTransient_h
#define genericSources_FitsTransient_h

#include "flux/Spectrum.h"

namespace IRB {
   class EblAtten;
}

/**
 * @class FitsTransient
 * @brief A flaring source whose spectral evolution is given by data in a 
 * FITS binary table.  The mean photon flux and flare start and stop times
 * are given as parameters.
 *
 * @author J. Chiang
 *
 * $Header$
 */

class FitsTransient : public Spectrum {

public:

   FitsTransient(const std::string & params);
   
   virtual ~FitsTransient();

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
      return "FitsTransient";
   }

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
   std::string m_fitsFile;
   float m_z;

   IRB::EblAtten * m_tau;
   double m_currentEnergy;
   std::vector<std::pair<double, double> > m_events;

   bool m_haveFirstEvent;
   std::vector<std::pair<double, double> >::const_iterator m_nextEvent;

   void createEvents();

   void readTable(const std::string & extname,
                  std::vector<double> & data,
                  const std::string & colname) const;

   void readSpectra(std::vector< std::vector<double> > & spectra) const;

   void computeIntegralDist(const std::vector<double> & times,
                            const std::vector<double> & energies,
                            const std::vector< std::vector<double> > & spectra,
                            std::vector< std::vector<double> > & integralDist,
                            std::vector<double> & lightCurve) const;

   void drawEvents(const std::vector<double> & times,
                   const std::vector<double> & energies,
                   const std::vector< std::vector<double> > & integralDist,
                   const std::vector<double> & lightCurve);

   std::pair<long, double> draw(const std::vector<double> & x,
                                const std::vector<double> & y) const;

};

#endif // genericSources_FitsTransient_h
