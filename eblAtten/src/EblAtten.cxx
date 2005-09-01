/**
 * @file EblAtten.cxx
 * @brief Function object wrapper to code in IRB_routines.cxx that
 * calculates EBL optical depth as a function of energy and redshift
 * for four different models.
 * @author J. Chiang
 *
 * $Header$
 */

#include <sstream>

#include "eblAtten/EblAtten.h"

namespace IRB {

std::map<EblModel, std::string> EblAtten::s_model_Ids;

float calcSS(float energy, float redshift);
float calcPB99(float energy, float redshift);
float calcKneiske(float energy, float redshift);

EblAtten::EblAtten(EblModel model) : m_model(model) {
   if (s_model_Ids.size() == 0) {
      s_model_Ids[Salamon_Stecker] = "Salamon & Stecker (1998)";
      s_model_Ids[Primack_Bullock99] = "Primack & Bullock (1999)";
      s_model_Ids[Kneiske] = "Kneiske et al. (2004)";
   }
   if (s_model_Ids.find(model) == s_model_Ids.end()) {
      std::ostringstream message;
      message << "Invalid model ID: " << model << "\n"
              << "Valid models are \n";
      std::map<EblModel, std::string>::iterator it;
      for (it = s_model_Ids.begin(); it != s_model_Ids.end(); ++it) {
         message << it->first << " : " << it->second << "\n";
      }
      throw std::runtime_error(message.str());
   }
}

float EblAtten::operator()(float energy, float redshift) const {
// Convert energy from MeV to TeV:
   energy /= 1e6;
   switch (m_model) {
   case Salamon_Stecker:
      return calcSS(energy, redshift);
   case Primack_Bullock99:
      return calcPB99(energy, redshift);
   case Kneiske:
      return calcKneiske(energy, redshift);
   }
   return 0;
}

} // namespace IRB
