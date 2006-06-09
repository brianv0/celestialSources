/** 
 * @file FileSpectrum.cxx
 * @brief Implementation of FileSpectrum
 *
 *  $Header$
 */

#include <cmath>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "CLHEP/Random/RandFlat.h"

#include "facilities/Util.h"

#include "flux/SpectrumFactory.h"

#include "genericSources/FileSpectrum.h"
#include "ConstParMap.h"
#include "Util.h"

namespace {
   double pl_integral(double e1, double e2, double f1, double f2) {
      double gamma(std::log(f2/f1)/std::log(e2/e1));
      double n0(f1/std::pow(e1, gamma));
      double my_value(0);
      if (gamma != -1) {
         my_value = n0/(gamma + 1.)*(std::pow(e2, gamma + 1.) - 
                                     std::pow(e1, gamma + 1.));
      } else {
         my_value = n0*std::log(e2/e1);
      }
      return my_value;
   }

   double pl_draw(double e1, double e2, double f1, double f2) {
      double gamma(std::log(f2/f1)/std::log(e2/e1));
      double e1gam(std::pow(e1, gamma));
      double e2gam(std::pow(e2, gamma));

      double xi(CLHEP::RandFlat::shoot());
      double ee(std::pow(xi*(e1gam - e2gam) + e2gam, 1./(gamma)));

      return ee;
   }
}

ISpectrumFactory & FileSpectrumFactory() {
   static SpectrumFactory<FileSpectrum> factory;
   return factory;
}

FileSpectrum::FileSpectrum(const std::string& params) 
   : m_emin(20), m_emax(2e5) {
   std::map<std::string, std::string> my_parmap;
   facilities::Util::keyValueTokenize(params, ",", my_parmap);
   
   genericSources::ConstParMap parmap(my_parmap);
   
   m_flux = parmap.value("flux");

   try {
      m_emin = parmap.value("emin");
   } catch (...) {
   }
   try {
      m_emax = parmap.value("emax");
   } catch (...) {
   }

   std::string infile = parmap["specFile"];
   double file_flux = read_file(infile);

// Set the flux to the file integral, if it is not set by the XML definition.
   if (m_flux == 0) {
      m_flux = file_flux;
   }
}

double FileSpectrum::flux() const {
   return m_flux;
}

double FileSpectrum::flux(double time) const {
   (void)(time);
   return m_flux;
}

float FileSpectrum::operator() (float xi) {
   size_t k(std::upper_bound(m_integralSpectrum.begin(),
                             m_integralSpectrum.end(), xi) 
            - m_integralSpectrum.begin());
   return ::pl_draw(m_energies.at(k-1), m_energies.at(k),
                    m_integralSpectrum.at(k-1), m_integralSpectrum.at(k));
}

std::string FileSpectrum::title() const {
    return "FileSpectrum";
}

const char * FileSpectrum::particleName() const {
    return m_particle_name.c_str();
}

double FileSpectrum::read_file(const std::string & infile) {
   genericSources::Util::file_ok(infile);
   std::vector<std::string> lines;
   genericSources::Util::readLines(infile, lines, "%");

   m_energies.clear();
   std::deque<double> dnde;

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
// Use a deque and push_front to get the spectra to go from high to low
// energies. See compute_integral_dist implementation.
      m_energies.push_front(std::atof(tokens.at(0).c_str()));
      dnde.push_front(std::atof(tokens.at(1).c_str()));
   }
   reset_ebounds(dnde);
   return compute_integral_dist(dnde);
}

void FileSpectrum::reset_ebounds(const std::deque<double> & dnde) {
   (void)(dnde);
// Here awaits a proper implementation by JCT.  For now, do nothing
// and just draw from the energy range given in the spectrum file.

//    if (m_emin < m_energies.front() && m_emax > m_energies.back()) {
// // Don't need to do anything, since bounds of m_energies are used 
// // automatically, so just print warning to the screen.
//       std::cout << "Resetting energy bounds to match those in the input "
//                 << "spectrum file." << std::endl;
//    } if (m_emin > m_energies.front()) {
//       std::cout << "Minimum energy in XML definition is greater than "
//                 << "than the minimum energy in the input spectrum file.\n"
//                 << "Resetting minimum energy bound to match the XML value."
//                 << std::endl;
//       std::vector<double>::iterator emin 
//          = std::upper_bound(m_energies.begin(), m_energies.end(), m_emin);
//       size_t k = emin - m_energies.begin();
//       std::vector<double> energies(m_energies.size() - k);
//       std::vector<double> new_dnde(energies.size());
//       energies.at(0) = m_emin;
//       new_dnde.at(0) = 
//       std::copy(emin, m_energies.end(), energies.begin() + 1);
}

double FileSpectrum::
compute_integral_dist(const std::deque<double> & dnde) {
// Form the integral distribution going from high energies to low
// energies since most spectra will be falling as a function of energy
// and we want rare, high energy events to be sampled near zero where
// sampling on a logrithmic grid is hard.
   m_integralSpectrum.clear();
   m_integralSpectrum.push_back(0);
   for (size_t k = 1; k < m_energies.size(); k++) {
      double dn(::pl_integral(m_energies.at(k-1), m_energies.at(k),
                              dnde.at(k-1), dnde.at(k)));
      m_integralSpectrum.push_back(m_integralSpectrum.back() + dn);
   }
   double total_flux(m_integralSpectrum.back());
   for (size_t k = 0; k < m_integralSpectrum.size(); k++) {
      m_integralSpectrum.at(k) /= total_flux;
   }

   std::ofstream outfile("integral_dist.dat");
   for (size_t k = 0; k < m_energies.size(); k++) {
      outfile << m_energies.at(k) << "  "
              << m_integralSpectrum.at(k) << std::endl;
   }
   outfile.close();

   return total_flux;
}
