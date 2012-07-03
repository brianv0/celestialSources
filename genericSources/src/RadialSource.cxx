/**
 * @file RadialSource.cxx
 * @brief A simple Spectrum subclass that exercises the flux package.
 * @author J. Chiang
 *
 * $Header$
 */

#include <cmath>

#include <algorithm>

#include "CLHEP/Random/RandFlat.h"

#include "astro/SkyDir.h"

#include "flux/SpectrumFactory.h"
#include "flux/EventSource.h"

#include "celestialSources/ConstParMap.h"

#include "genericSources/FileSpectrum.h"
#include "genericSources/RadialSource.h"

#include "Util.h"

ISpectrumFactory & RadialSourceFactory() {
   static SpectrumFactory<RadialSource> myFactory;
   return myFactory;
}

RadialSource::RadialSource(const std::string & paramString) 
   : m_spectrum(0) {

   celestialSources::ConstParMap pars(paramString);
   
// The following call requires the flux parameter to be set.  The actual
// value is held internally by the m_spectrum data member.
   pars.value("flux");

   m_spectrum = new FileSpectrum(paramString);

   fillRadialDist(pars["profileFile"]);

   double ra(pars.value("ra"));
   double dec(pars.value("dec"));
   astro::SkyDir srcDir(ra, dec);
   m_rot = 
      CLHEP::HepRotation().rotateZ(-ra*M_PI/180.).rotateY((dec-90.)*M_PI/180.);
   m_rotatedSrcVec = m_rot*srcDir.dir().unit();
}

float RadialSource::operator()(float xi) const {
   return (*m_spectrum)(xi);
}

double RadialSource::flux(double time) const {
   return m_spectrum->flux(time);
}

double RadialSource::solidAngle() const {
   return 1;
}

double RadialSource::interval(double time) {
   return m_spectrum->interval(time);
}

double RadialSource::energy(double time) {
   (void)(time);
   return (*m_spectrum)(CLHEP::RandFlat::shoot());
}

std::pair<double, double> RadialSource::dir(double energy) {
   (void)(energy);

   double phi(CLHEP::RandFlat::shoot()*2.*M_PI);
   double theta(drawOffsetAngle());

   CLHEP::Hep3Vector appDir 
      = CLHEP::HepRotation().rotateY(theta).rotateZ(phi)*m_rotatedSrcVec;

   astro::SkyDir myDir(m_rot.inverse()*appDir);
   
   return std::make_pair(myDir.l(), myDir.b());
}

void RadialSource::fillRadialDist(const std::string & infile) {
   genericSources::Util::file_ok(infile);
   std::vector<std::string> lines;
   genericSources::Util::readLines(infile, lines, "%#");
   m_thetas.clear();
   m_radialDist.clear();

   std::vector<double> dndth;

   std::vector<std::string>::const_iterator line = lines.begin();
   for ( ; line != lines.end(); ++line) {
      std::vector<std::string> tokens;
      facilities::Util::stringTokenize(*line, " \t", tokens);
      if (tokens.size() < 2) {
         std::ostringstream message;
         message << "RadialSource: poorly formatted column in input file: "
                 << infile;
         throw std::runtime_error(message.str());
      }
      m_thetas.push_back(std::atof(tokens.at(0).c_str())*M_PI/180.);
      dndth.push_back(std::atof(tokens.at(1).c_str()));
   }

   m_radialDist.push_back(0);
   for (size_t i(1); i < m_thetas.size(); i++) {
      double delta((dndth.at(i)*std::sin(m_thetas.at(i)) 
                    + dndth.at(i-1)*std::sin(m_thetas.at(i-1)))/2.
                   *(m_thetas.at(i) - m_thetas.at(i-1)));
      m_radialDist.push_back(m_radialDist.back() + delta);
   }
}

double RadialSource::drawOffsetAngle() const {
   double xi(CLHEP::RandFlat::shoot()*m_radialDist.back());
   std::vector<double>::const_iterator it
      = std::upper_bound(m_radialDist.begin(), m_radialDist.end(), xi);
   size_t indx = it - m_radialDist.begin() - 1;
   double theta = ( (xi - m_radialDist.at(indx))
                    /(m_radialDist.at(indx+1) - m_radialDist.at(indx))
                    *(m_thetas.at(indx+1) - m_thetas.at(indx)) 
                    + m_thetas.at(indx) );
   return theta;
}
