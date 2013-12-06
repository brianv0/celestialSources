/**
 * @file AsciiTableModel.h
 * @brief Interface to ascii tabulation of EBL optical depth vs 
 * (redshift, energy) using the same format as the Gilmore 2012 optical 
 * depth tables.
 * 
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef IRB_AsciiTableModel_h
#define IRB_AsciiTableModel_h

#include <string>
#include <vector>

namespace IRB {

class AsciiTableModel {

public:

   AsciiTableModel(const std::string & infile);

   float value(float energy, float redshift) const;

private:

   std::string m_infile;
   std::vector<float> m_energies;
   std::vector<float> m_redshifts;
   std::vector< std::vector<float> > m_tau_array;

   void read_ascii_table();

   float tau_of_e(float redshift, size_t e_index, size_t z_index) const;

};

} // namespace IRB

#endif // IRB_AsciiTableModel_h
