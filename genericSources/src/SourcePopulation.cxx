/**
 * @file SourcePopulation.cxx
 * @brief Generate events from a population of steady point sources.
 * @author J. Chiang
 *
 * $Header$
 */

#include <cmath>
#include <cstdlib>

#include <algorithm>
#include <fstream>
#include <iostream>

#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoisson.h"

#include "facilities/Util.h"

#include "astro/SkyDir.h"

#include "flux/SpectrumFactory.h"

#include "eblAtten/EblAtten.h"

#include "genericSources/SourcePopulation.h"

namespace {
   void cleanLine(std::string & line) {
      char CR[1];
      CR[0] = 0x0d;
      if (line.find(CR) != std::string::npos) {
         std::vector<std::string> tokens;
         facilities::Util::stringTokenize(line, CR, tokens);
         line = tokens.front();
      }
   }
   
   void readLines(std::string inputFile, 
                  std::vector<std::string> & lines,
                  const std::string & skip,
                  bool cleanLines) {
      facilities::Util::expandEnvVar(&inputFile);
      std::ifstream file(inputFile.c_str());
      lines.clear();
      std::string line;
      while (std::getline(file, line, '\n')) {
         if (line != "" && line != " "             //skip (most) blank lines 
             && line.find_first_of(skip) != 0) {   //and commented lines
            if (cleanLines) {
               cleanLine(line);
            }
            lines.push_back(line);
         }
      }
   }
}

ISpectrumFactory & SourcePopulationFactory() {
   static SpectrumFactory<SourcePopulation> myFactory;
   return myFactory;
}

SourcePopulation::SourcePopulation(const std::string & params) {
   std::vector<std::string> pars;
   facilities::Util::stringTokenize(params, ",", pars);
   std::string inputFile(pars.at(0));
   facilities::Util::expandEnvVar(&inputFile);
   readSourceFile(inputFile);
   setFlux(m_cumulativeFlux.back());
}

void SourcePopulation::readSourceFile(const std::string & input_file) {
   std::vector<std::string> lines;
   ::readLines(input_file, lines, "#", true);

   m_sources.clear();
   m_sources.reserve(lines.size());

   m_cumulativeFlux.clear();
   m_cumulativeFlux.reserve(lines.size());

   std::vector<std::string>::const_iterator line = lines.begin();
   for ( ; line != lines.end(); ++line) {
      std::vector<std::string> tokens;
      facilities::Util::stringTokenize(*line, ", \n", tokens);
      double ra = std::atof(tokens.at(0).c_str());
      double dec = std::atof(tokens.at(1).c_str());
      double flux = std::atof(tokens.at(2).c_str());
      double gamma = std::atof(tokens.at(3).c_str());
      double gamma2 = std::atof(tokens.at(4).c_str());
      double ebreak = std::atof(tokens.at(5).c_str());
      double emin = std::atof(tokens.at(6).c_str());
      double emax = std::atof(tokens.at(7).c_str());
      m_sources.push_back(PointSource(astro::SkyDir(ra, dec), flux, gamma,
                                      gamma2, ebreak, emin, emax));
      if (line == lines.begin()) {
         m_cumulativeFlux.push_back(flux);
      } else {
         m_cumulativeFlux.push_back(m_cumulativeFlux.back() + flux);
      }
   }
}

double SourcePopulation::interval(double time) {
   (void)(time);
   double xi(RandFlat::shoot()*m_flux);
   std::vector<double>::const_iterator it = 
      std::upper_bound(m_cumulativeFlux.begin(), m_cumulativeFlux.end(), xi);
   size_t indx = it - m_cumulativeFlux.begin();
   m_currentEnergy = m_sources.at(indx).energy();
   m_l = m_sources.at(indx).dir().l();
   m_b = m_sources.at(indx).dir().b();
   return -1;
}

SourcePopulation::
PointSource::PointSource(const astro::SkyDir & dir, double flux,
                         double gamma, double gamma2, double ebreak,
                         double emin, double emax, const IRB::EblAtten * tau) 
   : m_dir(dir), m_flux(flux), m_gamma(gamma), m_gamma2(gamma2),
     m_ebreak(ebreak), m_emin(emin), m_emax(emax), m_tau(tau) {

   if (gamma == 1) {
      m_part1 = std::log(ebreak/emin);
   } else {
      m_part1 = (std::pow(ebreak, 1 - gamma) - std::pow(emin, 1 - gamma))
         /(1 - gamma)/std::pow(ebreak, 1 - gamma);
   }
   if (gamma2 == 1) {
      m_part2 = std::log(emax/ebreak);
   } else {
      m_part2 = (std::pow(emax, 1 - gamma2) - std::pow(ebreak, 1 - gamma2))
         /(1 - gamma2)/std::pow(ebreak, 1 - gamma2);
   }
   m_frac = m_part1/(m_part1 + m_part2);
}

SourcePopulation::
PointSource::~PointSource() {
   delete m_tau;
}

double SourcePopulation::
PointSource::energy() const {
   double xi(RandFlat::shoot());
   if (xi < m_frac) {
      xi /= m_frac;
      double aa(1. - m_gamma);
      double bb(std::pow(m_ebreak, 1. - m_gamma));
      return std::pow(xi*aa*bb*m_part1 + std::pow(m_emin, aa), 1./aa);
   }
   xi = (xi - m_frac)/(1. - m_frac);
   double aa(1. - m_gamma2);
   double bb(std::pow(m_ebreak, 1. - m_gamma2));
   return std::pow(xi*aa*bb*m_part2 + std::pow(m_ebreak, aa), 1./aa);
}
