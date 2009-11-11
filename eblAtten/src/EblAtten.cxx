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
float calcFranceschini(float energy, float redshift);
float calcFinke(float energy, float redshift);
float calcGilmore(float energy, float redshift);
float calcStecker05_FE(float energy, float redshift);
float calcSalamonStecker(float energy, float redshift);
float calcGeneric(float energy, float redshift);

EblAtten::EblAtten(EblModel model) : m_model(model) {
   if (s_model_Ids.size() == 0) {
      s_model_Ids[Kneiske] = "Kneiske et al - Best Fit (2004)";
      s_model_Ids[Primack05] = "Primack et al (2005)";
      s_model_Ids[Kneiske_HighUV] = "Kneiske et al - High UV (2004)";
      s_model_Ids[Stecker05] = "Stecker et al (2006)";
      s_model_Ids[Franceschini] = "Franceschini (2008)";
      s_model_Ids[Finke] = "Finke et al. (2009)";
      s_model_Ids[Gilmore] = "Gilmore et al. (2008)";
      s_model_Ids[Stecker05_FE] = "Stecker et al (2006) - Fast Evolution";
      s_model_Ids[SalamonStecker] = "Salamon & Stecker (1998) - with metallicity correction";
      s_model_Ids[Generic] = "Generic representation of tau(E,z) from Justin Finke"; 	  
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
   case Franceschini:
      return calcFranceschini(energy, redshift);
   case Finke:
      return calcFinke(energy, redshift);
   case Gilmore:
      return calcGilmore(energy, redshift);
   case Stecker05_FE:
      return calcStecker05_FE(energy, redshift);
   case SalamonStecker:
      return calcSalamonStecker(energy, redshift);
   case Generic:
      return calcGeneric(energy, redshift);
   }
   return 0;
}

} // namespace IRB
