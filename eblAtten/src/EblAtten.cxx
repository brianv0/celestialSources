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
float calcGilmore12_fixed(float energy, float redshift);
float calcGilmore12_fiducial(float energy, float redshift);
float calcInoue13(float energy, float redshift);
float calcDominguez11(float energy, float redshift);
float calcScully14_highOp(float energy, float redshift);
float calcScully14_lowOp(float energy, float redshift);
float calcKneiskeDole10(float energy, float redshift);
float calcKneiskeDole10_CMB(float energy, float redshift);

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
      s_model_Ids[Gilmore12_fixed] = "Gilmore et al (2012) - WMAP5+Fixed";
      s_model_Ids[Gilmore12_fiducial] = "Gilmore et al (2012) - Evolving Dust";
      s_model_Ids[Inoue13] = "Inoue et al (2013)";
      s_model_Ids[Dominguez11] = "Dominguez et al (2011)";
      s_model_Ids[Scully14_highOp] = "Scully et al (2014) - High Opacity";
      s_model_Ids[Scully14_lowOp] = "Scully et al (2014) - Low Opacity";
      s_model_Ids[KneiskeDole10] = "Kneiske & Dole (2010)";
      s_model_Ids[KneiskeDole10_CMB] = "Kneiske & Dole (2010) - with CMB";
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
   case Gilmore12_fixed:
      return calcGilmore12_fixed(energy, redshift);
   case Gilmore12_fiducial:
      return calcGilmore12_fiducial(energy, redshift);
   case Inoue13:
      return calcInoue13(energy, redshift);
   case Dominguez11:
      return calcDominguez11(energy, redshift);
   case Scully14_highOp:
      return calcScully14_highOp(energy, redshift);
   case Scully14_lowOp:
      return calcScully14_lowOp(energy, redshift);
   case KneiskeDole10:
      return calcKneiskeDole10(energy, redshift);
   case KneiskeDole10_CMB:
      return calcKneiskeDole10_CMB(energy, redshift);
   }
   return 0;
}

} // namespace IRB
