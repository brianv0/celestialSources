/** 
 * @file FileSpectrum.cxx
 * @brief Implementation of FileSpectrum
 *
 *  $Header$
 */

#include <cmath>

#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "facilities/Util.h"

#include "flux/SpectrumFactory.h"

#include "genericSources/FileSpectrum.h"

ISpectrumFactory & FileSpectrumFactory() {
   static SpectrumFactory<FileSpectrum> factory;
   return factory;
}

FileSpectrum::FileSpectrum(const std::string& params) {
   std::vector<EFpair_t> temp_vector;
// retrieve and parse the parameter string
   std::map<std::string, std::string> parmap;
   facilities::Util::keyValueTokenize(params, ",", parmap);

   m_flux = std::atof(parmap["flux"].c_str());

// open the spectrum file
   std::string infile = parmap["specFile"];
   double file_flux = read_file(infile);

// set the flux to the file integral, if not set by XML
   if (m_flux == 0) {
      m_flux = file_flux;
   }
}

double FileSpectrum::flux() const {
   return m_flux;
}

double FileSpectrum::flux(double time) const {
   return m_flux;
}

float FileSpectrum::operator() (float xi) {
   std::vector<EFpair_t>::const_iterator it
      = std::upper_bound(m_integralSpectrum.begin(), m_integralSpectrum.end(),
                         std::make
}

std::string FileSpectrum::title() const {
    return "FileSpectrum";
}

const char * FileSpectrum::particleName() const {
    return m_particle_name.c_str();
}

double FileSpectrum::read_file(const std::string & infile) {
   Util::file_ok(infile);
   std::vector<std::string> lines;
   Util::readLines(infile, lines, "%");

   std::vector<double> energies;
   std::vector<double> dNdE;

   std::vector<std::string>::const_iterator line = lines.begin();
   for ( ; line != lines.end(); ++line) {
      std::vector<std::string> tokens;
      facilities::Util::stringTokenize(*line, " \t", tokens);
      if (tokens.size() < 2) {
         std::ostringstream message;
         message << "FileSpectrum: poorly formatted column in input file: "
                 << infile;
         throw std::runtime_error(message.str());
      }
      energies.push_back(std::atof(tokens.at(0).c_str()));
      fluxes.push_back(std::atof(tokens.at(1).c_str()));
   }
   return compute_integral_dist(energies, fluxes);
}

double FileSpectrum::
compute_integral_dist(const std::vector<double> & energies,
                      const std::vector<double> & fluxes) {
   m_integralSpectrum.clear();
   m_integralSpectrum.reserve(energies.size());
   m_integralSpectrum.push_back(std::make_pair(energies.front(), 0));
   for (size_t k = 1; k < energies.size(); k++) {
      double dflux = ( (fluxes.at(k) + fluxes.at(k-1))/2.
                       *(energies.at(k) - energies.at(k-1)) );
      double integral(m_integralSpectrum.back(k).second + dflux);
      m_integralSpectrum.push_back(std::make_pair(energies.at(k), integral));
   }
   double total_flux(m_integralSpectrum.back().second);
   for (size_t k = 0; k < m_integralSpectrum.size(); k++) {
      m_integralSpectrum.at(k).second /= total_flux;
   }
   return total_flux;
}
