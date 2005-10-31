/**
 * @file MapCube.h
 * @brief A source class for the flux package that uses FITS images as
 * templates for the photon distribution on the sky.
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef mySpectrum_MapCube_h
#define mySpectrum_MapCube_h

#include <vector>

#include "genericSources/MapSource.h"

#include "fitsio.h"

/**
 * @class MapCube
 *
 * @brief A source class for the flux package that uses FITS images as
 * templates for the photon distribution on the sky.  The image data
 * should comprise data cube with two spatial dimensions and an energy
 * dimension.
 *
 * @author J. Chiang
 *
 * $Header$
 */

class MapCube : public MapSource {

public:

   /// This constructor is required and used in FluxSource for
   /// "SpectrumClass" sources.
   MapCube(const std::string &params);

   virtual ~MapCube() {}

   /// @return Particle energy in MeV.
   /// @param xi Uniform random deviate on the unit interval.
   virtual float operator() (float xi) const;

   /// @return Photon energy (MeV).
   virtual double energy(double time);

   /// @return Photon direction in (l, b).
   /// @param energy This is the energy returned from a previous call
   /// to the energy(time) method (see the flux/ISpectrum
   /// declaration).  For diffuse sources with spatially varying
   /// spectra, this logic is *backwards*.  One would want to find the
   /// incident direction first, then find the energy.
   virtual std::pair<double, double> dir(double energy) {
      (void)(energy);
      return m_currentDir;
   }

   /// @return Particle type, "gamma".
   virtual const char * particleName() const {return "gamma";}

   /// @return Title describing the spectrum.
   virtual std::string title() const {return "MapCube";}

private:

   /// This is mutable since operator()(float) is const (inherited
   /// from the flux/Spectrum base class).
   mutable std::pair<double, double> m_currentDir;

   std::vector<std::vector<std::pair<double,double> > > m_spectra;
   std::vector<double> m_energies;

   double mapValue(unsigned int i, unsigned int j, unsigned int k);

   void readEnergyVector(const std::string & fitsFile);
   
   void readColumn(fitsfile * fptr, const std::string & colname,
                   std::vector<double> & coldata) const;

   void makeCumulativeSpectra();

   double powerLawIntegral(double x1, double x2, double y1, double y2,
                           double & gamma) const;
                           
   double drawEnergy(const std::vector<std::pair<double, double> > & spectrum)
      const;

   static bool cmpPair(const std::pair<double, double> & x, 
                       const std::pair<double, double> & y); 
};

#endif // mySpectrum_MapCube_h
