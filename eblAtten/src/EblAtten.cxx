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

float calcKneiske(float energy, float redshift);
float calcPrimack05(float energy, float redshift);
float calcKneiske_HighUV(float energy, float redshift);
float calcStecker05(float energy, float redshift);

EblAtten::EblAtten(EblModel model) : m_model(model) {
   if (s_model_Ids.size() == 0) {
      s_model_Ids[Kneiske] = "Kneiske et al - Best Fit (2004)";
      s_model_Ids[Primack05] = "Primack et al (2005)";
      s_model_Ids[Kneiske_HighUV] = "Kneiske et al - High UV (2004)";
      s_model_Ids[Stecker05] = "Stecker et al (2005)";
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
// Convert energy from MeV to GeV:
   energy /= 1e3;
   switch (m_model) {
   case Kneiske:
      return calcKneiske(energy, redshift);
   case Primack05:
      return calcPrimack05(energy, redshift);
   case Kneiske_HighUV:
      return calcKneiske_HighUV(energy, redshift);
   case Stecker05:
      return calcStecker05(energy, redshift);
   }
   return 0;
}

} // namespace IRB
