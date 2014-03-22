/**
 * @file AsciiTableModel.cxx
 * @brief Interface to ascii tabulation of EBL optical depth vs 
 * (redshift, energy) using the same format as the Gilmore 2012 optical 
 * depth tables.
 * 
 * @author J. Chiang
 *
 * $Header$
 */

#include <cmath>

#include <algorithm>
#include <stdexcept>

#include "facilities/Util.h"
#include "st_facilities/Util.h"
#include "AsciiTableModel.h"

namespace IRB {

AsciiTableModel::AsciiTableModel(const std::string & infile) 
   : m_infile(infile) {
   read_ascii_table();
}

float AsciiTableModel::value(float energy, float redshift) const {
   if (redshift < m_redshifts.front() || energy < m_energies.front()) {
      return 0;
   }

   size_t e_index = std::lower_bound(m_energies.begin(), m_energies.end(),
                                     energy) - m_energies.begin() - 1;


   if (e_index + 1  > m_energies.size() - 1) {
      throw std::runtime_error("Selected energy outside range of "
                               + m_infile);
   }
   size_t z_index = std::lower_bound(m_redshifts.begin(), m_redshifts.end(),
                                     redshift) - m_redshifts.begin() -1 ;

   if (z_index + 1 > m_redshifts.size() - 1) {
      throw std::runtime_error("Selected redshift outside range of "
                               + m_infile);
   }
   float tau1 = tau_of_e(redshift, e_index, z_index);
   float tau2 = tau_of_e(redshift, e_index+1, z_index);

   // Interpolate in log-log space.
   float emin(m_energies[e_index]);
   float emax(m_energies[e_index + 1]);
   if (tau1 == 0 || tau2 == 0) {
      return 0;
   }
   float tau = ( tau1*std::exp(std::log(energy/emin)/std::log(emax/emin)
                               *std::log(tau2/tau1)) );
   return tau;
}

void AsciiTableModel::read_ascii_table() {
   std::vector<std::string> lines;
   std::string skip;
   bool cleanLines;
   st_facilities::Util::readLines(m_infile, lines, skip="", cleanLines=true);

   // Read the redshift values from the first line.
   std::vector<std::string> tokens;
   facilities::Util::stringTokenize(lines[0].substr(16), ",", tokens);
   m_redshifts.clear();
   for (size_t i(0); i < tokens.size(); i++) {
      m_redshifts.push_back(std::atof(tokens[i].c_str()));
   }
   // Read in the energies and tau values.
   m_energies.clear();
   m_tau_array.clear();
   for (size_t i(1); i < lines.size(); i++) {
      facilities::Util::stringTokenize(lines[i], "\t ", tokens);
      if (tokens.size() != m_redshifts.size() + 1) {
         continue;
      }
      m_energies.push_back(std::atof(tokens[0].c_str()));
      std::vector<float> tau_values;
      for (size_t k(1); k < tokens.size(); k++) {
         tau_values.push_back(std::atof(tokens[k].c_str()));
      }
      m_tau_array.push_back(tau_values);
   }
}

float AsciiTableModel::tau_of_e(float redshift, size_t e_index,
                                size_t z_index) const {
   float tau = ( (redshift - m_redshifts[z_index])
                 /(m_redshifts[z_index+1] - m_redshifts[z_index])
                 *(m_tau_array[e_index][z_index+1] - m_tau_array[e_index][z_index])
                 + m_tau_array[e_index][z_index] );
   return tau;
}

} // namespace IRB
