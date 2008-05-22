/**
 * @file FitsTransient.cxx
 * @brief A flaring source whose spectral evolution is given by a FITS
 * binary table.
 * @author J. Chiang
 *
 * $Header$
 */

#include <algorithm>
#include <iostream>
#include <map>
#include <utility>
#include <vector>

#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoisson.h"

#include "facilities/Util.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "flux/EventSource.h"
#include "flux/SpectrumFactory.h"

#include "eblAtten/EblAtten.h"

#include "genericSources/FitsTransient.h"

#include "ConstParMap.h"
#include "Util.h"

namespace {
   bool compareEventTime(const std::pair<double, double> & x,
                         const std::pair<double, double> & y) {
      return x.first < y.first;
   }
}

ISpectrumFactory & FitsTransientFactory() {
   static SpectrumFactory<FitsTransient> myFactory;
   return myFactory;
}

FitsTransient::FitsTransient(const std::string & paramString) 
   : m_z(0), m_tau(0), m_haveFirstEvent(false) {
   if (paramString.find("=") == std::string::npos) {
      std::vector<std::string> params;
      facilities::Util::stringTokenize(paramString, ", ", params);
      
      m_flux = std::atof(params[0].c_str());
      m_tstart = std::atof(params[1].c_str());
      m_tstop = std::atof(params[2].c_str());
      m_fitsFile = params[3];
      if (params.size() > 4) {
         m_z = std::atof(params[4].c_str());
      }
      if (params.size() > 5) {
         IRB::EblModel eblModel = 
            static_cast<IRB::EblModel>(std::atoi(params[5].c_str()));
         try {
            m_tau = new IRB::EblAtten(eblModel);
         } catch (std::exception & eObj) {
            std::cerr << eObj.what() << "\n"
                      << "Using default, IRB::Kneiske" << std::endl;
            m_tau = new IRB::EblAtten(IRB::Kneiske);
         }
      } else {
         m_tau = new IRB::EblAtten(IRB::Kneiske);
      }
   } else {
      genericSources::ConstParMap parmap(paramString);
      m_flux = parmap.value("flux");
      m_tstart = parmap.value("tstart");
      m_tstop = parmap.value("tstop");
      m_fitsFile = parmap["fitsFile"];
      try {
         m_z = parmap.value("z");
         m_tau = new IRB::EblAtten(IRB::Kneiske);
      } catch (...) {
      }
      try {
         IRB::EblModel eblModel = 
            static_cast<IRB::EblModel>(std::atoi(parmap["eblModel"].c_str()));
         try {
            m_tau = new IRB::EblAtten(eblModel);
         } catch (std::exception & eObj) {
            std::cerr << eObj.what() << "\n"
                      << "Using default, IRB::Kneiske" << std::endl;
            m_tau = new IRB::EblAtten(IRB::Kneiske);
         }
      } catch (...) {
         m_tau = new IRB::EblAtten(IRB::Kneiske);
      }
   }

   createEvents();
}

FitsTransient::~FitsTransient() {
   delete m_tau;
}

double FitsTransient::interval(double time) {
   time -= Spectrum::startTime();
   if (!m_haveFirstEvent) {
      std::vector<std::pair<double, double> >::const_iterator event =
         std::upper_bound(m_events.begin(), m_events.end(), 
                          std::make_pair(time, 0), compareEventTime);
      m_nextEvent = event;
      m_haveFirstEvent = true;
   }
   if (m_nextEvent != m_events.end()) {
      double dt(m_nextEvent->first - time);
      m_currentEnergy = m_nextEvent->second;
      ++m_nextEvent;
      return dt;
   }
// There should be a better way to turn off a source than this:
   return 3.155e8;
}

void FitsTransient::createEvents() {
   facilities::Util::expandEnvVar(&m_fitsFile);
   genericSources::Util::file_ok(m_fitsFile);

   std::vector<double> energies;
   std::vector<double> times;
   std::vector< std::vector<double> > spectra;
   readTable("ENERGIES", energies, "Energy");
   readTable("TIMES", times, "Time");
   double scale_factor((m_tstop - m_tstart)/(times.back() - times.front()));
   double t0(times.front());
   for (size_t i = 0; i < times.size(); i++) {
      times.at(i) = (times.at(i) - t0)*scale_factor + m_tstart;
   }
   readSpectra(spectra);
   std::vector< std::vector<double> > integralDist;
   std::vector<double> lightCurve;
   computeIntegralDist(times, energies, spectra, integralDist, lightCurve);
   drawEvents(times, energies, integralDist, lightCurve);
   std::stable_sort(m_events.begin(), m_events.end(), compareEventTime);
//    std::cout << "nevents = " << m_events.size() << std::endl;
//    for (size_t i(0); i < m_events.size(); i++) {
//       std::cout << m_events.at(i).first << "  "
//                 << m_events.at(i).second << std::endl;
//    }
}

void FitsTransient::
drawEvents(const std::vector<double> & times, 
           const std::vector<double> & energies,
           const std::vector< std::vector<double> > & integralDist,
           const std::vector<double> & lightCurve) {
   double npred = m_flux*EventSource::totalArea()*(m_tstop - m_tstart);
   long nevts = CLHEP::RandPoisson::shoot(npred);
   for (long i = 0; i < nevts; i++) {
      std::pair<long, double> arrTime = draw(times, lightCurve);
      std::pair<long, double> energy = draw(energies,
                                            integralDist.at(arrTime.first));
      m_events.push_back(std::make_pair(arrTime.second, energy.second));
   }
}

std::pair<long, double> FitsTransient::
draw(const std::vector<double> & x,
     const std::vector<double> & integralDist) const {
   double xi = CLHEP::RandFlat::shoot()*integralDist.back();
   std::vector<double>::const_iterator it = 
      std::upper_bound(integralDist.begin(), integralDist.end(), xi);
   long indx = it - integralDist.begin() - 1;
   double value = (x.at(indx+1) - x.at(indx))*(xi - integralDist.at(indx))/
      (integralDist.at(indx+1) - integralDist.at(indx)) + x.at(indx);
   return std::make_pair(indx, value);
}

void FitsTransient::
computeIntegralDist(const std::vector<double> & times,
                    const std::vector<double> & energies,
                    const std::vector< std::vector<double> > & spectra,
                    std::vector< std::vector<double> > & integralDist,
                    std::vector<double> & lightCurve) const {
   integralDist.resize(spectra.size());
   for (size_t i = 0; i < spectra.size(); i++) {
      integralDist.at(i).clear();
      integralDist.at(i).reserve(spectra.at(i).size());
      integralDist.at(i).push_back(0);
      for (size_t k = 1; k < spectra.at(i).size(); k++) {
         double counts = (spectra.at(i).at(k) + spectra.at(i).at(k-1))/2.
            *(energies.at(k) - energies.at(k-1));
         integralDist.at(i).push_back(integralDist.at(i).at(k-1) + counts);
      }
   }
   lightCurve.clear();
   lightCurve.push_back(0);
   for (size_t i = 1; i < spectra.size(); i++) {
      double cnts = (integralDist.at(i).back()+integralDist.at(i-1).back())/2.
         *(times.at(i) - times.at(i-1));
      lightCurve.push_back(lightCurve.at(i-1) + cnts);
   }
}

void FitsTransient::
readSpectra(std::vector< std::vector<double> > & spectra) const {
   std::auto_ptr<const tip::Table> 
      table(tip::IFileSvc::instance().readTable(m_fitsFile, "SPECTRA"));
   int nrows(table->getNumRecords());

   spectra.resize(nrows);

   tip::Table::ConstIterator row(table->begin());
   tip::ConstTableRecord & record(*row);

   std::vector< std::vector<double> >::iterator spectrum;
   for (spectrum = spectra.begin(); spectrum != spectra.end();
        ++spectrum, ++row) {
      record["Photon Spectrum"].get(*spectrum);
   }
}

void FitsTransient::readTable(const std::string & ext, 
                              std::vector<double> & data,
                              const std::string & colname) const {
   std::auto_ptr<const tip::Table> 
      table(tip::IFileSvc::instance().readTable(m_fitsFile, ext));
   int nrows(table->getNumRecords());

   data.resize(nrows);

   tip::Table::ConstIterator row(table->begin());
   tip::ConstTableRecord & record(*row);

   std::vector<double>::iterator datum;
   for (datum = data.begin(); datum != data.end(); ++datum, ++row) {
      record[colname].get(*datum);
   }
}
