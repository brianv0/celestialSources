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

#include "flux/EventSource.h"
#include "flux/SpectrumFactory.h"

#include "eblAtten/EblAtten.h"

#include "genericSources/SourcePopulation.h"

#include "ConstParMap.h"
#include "Util.h"
#include "GaussianQuadrature.h"

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

IRB::EblAtten * SourcePopulation::PointSource::s_tau(0);

ISpectrumFactory & SourcePopulationFactory() {
   static SpectrumFactory<SourcePopulation> myFactory;
   return myFactory;
}

SourcePopulation::SourcePopulation(const std::string & params) 
   : m_tau(0), m_idOffset(100000), m_l(0), m_b(0), m_name("") {
   if (params.find("=") == std::string::npos) {
      std::vector<std::string> pars;
      facilities::Util::stringTokenize(params, ",", pars);
      readSourceFile(pars.at(0));
      if (pars.size() > 1) {
         setEblAtten(pars.at(1));
      }
      if (pars.size() > 2) {
         m_idOffset = std::atoi(pars.at(2).c_str());
      } 
   } else {
      genericSources::ConstParMap pars(params);
      readSourceFile(pars["sourceFile"]);
      try {
         setEblAtten(pars["eblModel"]);
      } catch(std::runtime_error & ) {
      }
      try {
         m_idOffset = static_cast<int>(pars.value("idOffset"));
      } catch(std::runtime_error & ) {
      }
   }
   m_flux = m_cumulativeFlux.back();
}

SourcePopulation::~SourcePopulation() {
  delete m_tau;
}

void SourcePopulation::setEblAtten(const std::string & ebl_par) {
   IRB::EblModel eblModel =
      static_cast<IRB::EblModel>(std::atoi(ebl_par.c_str()));
   m_tau = new IRB::EblAtten(eblModel);
   PointSource::setEblAtten(m_tau);
}

void SourcePopulation::readSourceFile(std::string input_file) {
   facilities::Util::expandEnvVar(&input_file);
   genericSources::Util::file_ok(input_file);
   std::vector<std::string> lines;
   ::readLines(input_file, lines, "#", true);

   m_sources.clear();
   m_sources.reserve(lines.size());

   m_cumulativeFlux.clear();
   m_cumulativeFlux.reserve(lines.size());

   std::vector<std::string>::const_iterator line = lines.begin();
   for ( ; line != lines.end(); ++line) {
      m_sources.push_back(PointSource(*line));
      if (line == lines.begin()) {
         m_cumulativeFlux.push_back(m_sources.back().flux());
      } else {
         m_cumulativeFlux.push_back(m_cumulativeFlux.back() 
                                    + m_sources.back().flux());
      }
   }
}

float SourcePopulation::operator()(float xi) {
   xi *= m_flux;
   std::vector<double>::const_iterator it = 
      std::upper_bound(m_cumulativeFlux.begin(), m_cumulativeFlux.end(), xi);
   size_t indx = it - m_cumulativeFlux.begin();
   m_currentEnergy = m_sources.at(indx).energy();
   m_l = m_sources.at(indx).dir().l();
   m_b = m_sources.at(indx).dir().b();
   setIdentifier(indx + m_idOffset);
   m_name = m_sources.at(indx).name();
   return m_currentEnergy;
}

double SourcePopulation::energy(double time) {
   (void)(time);
   double xi = CLHEP::RandFlat::shoot();
   return this->operator()(xi);
}

double SourcePopulation::interval(double time) {
   (void)(time);
   double my_interval = 
      -std::log(1. - CLHEP::RandFlat::shoot())/m_flux/EventSource::totalArea();
   if (my_interval == 0) {
      throw std::runtime_error("SourcePopulation::interval:\n"
                               "zero length interval generated.");
   }
   return my_interval;
}

std::string SourcePopulation::name() const {
   return m_name;
}

//SourcePopulation::PointSource * SourcePopulation::PointSource::Self::s_self(0);

SourcePopulation::
PointSource::PointSource(const astro::SkyDir & dir, double flux,
                         double gamma, double gamma2, double ebreak,
                         double emin, double emax, double zz) 
   : m_dir(dir), m_flux(flux), m_gamma(gamma), m_gamma2(gamma2),
     m_ebreak(ebreak), m_emin(emin), m_emax(emax), m_z(zz) {
   setPowerLaw();
}

SourcePopulation::
PointSource::PointSource(const std::string & line) : m_z(0) {
   std::vector<std::string> tokens;
   facilities::Util::stringTokenize(line, ", \n", tokens);
   m_name = tokens.at(0);
   double ra = std::atof(tokens.at(1).c_str());
   double dec = std::atof(tokens.at(2).c_str());
   m_dir = astro::SkyDir(ra, dec);
   m_flux = std::atof(tokens.at(3).c_str());
   m_gamma = std::atof(tokens.at(4).c_str());
   m_gamma2 = std::atof(tokens.at(5).c_str());
   m_ebreak = std::atof(tokens.at(6).c_str());
   m_emin = std::atof(tokens.at(7).c_str());
   m_emax = std::atof(tokens.at(8).c_str());
   if (tokens.size() == 10) {
      m_z = std::atof(tokens.at(9).c_str());
   }
   setPowerLaw();
}

void 
SourcePopulation::
PointSource::setPowerLaw() {
   m_part1 = integral(m_emin, m_ebreak);
   m_part2 = integral(m_ebreak, m_emax);
   m_frac = m_part1/(m_part1 + m_part2);
}

double SourcePopulation::
PointSource::energy() const {
   double my_energy;
   double xi(CLHEP::RandFlat::shoot());
   if (xi < m_frac) {
      do {
         xi = CLHEP::RandFlat::shoot();
         double aa(1. - m_gamma);
         double bb(std::pow(m_ebreak, aa) - std::pow(m_emin, aa));
         my_energy = std::pow(xi*bb + std::pow(m_emin, aa), 1./aa);
      } while (CLHEP::RandFlat::shoot() > attenuation(my_energy));
   } else {
      do {
         xi = CLHEP::RandFlat::shoot();
         double aa(1. - m_gamma2);
         double bb(std::pow(m_emax, aa) - std::pow(m_ebreak, aa));
         my_energy = std::pow(xi*bb + std::pow(m_ebreak, aa), 1./aa);
      } while (CLHEP::RandFlat::shoot() > attenuation(my_energy));
   }
   return my_energy;
}

double SourcePopulation::
PointSource::integral(double emin, double emax) {
   double err(1e-5);
   double result(0);
   int ierr;
   DndeIntegrand dnde(*this);
   result = genericSources::GaussianQuadrature::dgaus8(dnde, emin, emax, 
                                                       err, ierr);
   return result;
}

double SourcePopulation::
PointSource::dnde(double energy) const {
   if (energy < m_emin || energy > m_emax) {
      return 0;
   }
   double gamma;
   if (energy < m_ebreak) {
      gamma = m_gamma;
   } else {
      gamma = m_gamma2;
   }
   return attenuation(energy)*std::pow(energy/m_ebreak, -gamma);
}

double SourcePopulation::
PointSource::attenuation(double energy) const {
   double atten(1.);
   if (s_tau != 0) {
      atten = std::exp(-s_tau->operator()(energy, m_z));
   }
   return atten;
}

