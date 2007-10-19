/**
 * @file EblAtten.cxx
 * @brief Function object wrapper to Hays/McEnery code (in IRB_routines.cxx) 
 * that calculates EBL optical depth as a function of energy and redshift
 * for four different models.
 * @author J. Chiang
 *
 * $Header$
 */

#include "genericSources/EblAtten.h"

namespace IRB {

float calcSJ1(float energy, float redshift);
float calcSJ2(float energy, float redshift);
float calcPrimack(float energy, float redshift);
float calcPrimack04(float energy, float redshift);

EblAtten::EblAtten(EblModel model) : m_model(model) {
   if (m_model < SdJbase || m_model > Primack04) {
      throw std::runtime_error("Invalid model.");
   }
}

float EblAtten::operator()(float energy, float redshift) const {
// Convert energy from MeV to TeV:
   energy /= 1e6;
   switch (m_model) {
   case SdJbase:
      return calcSJ1(energy, redshift);
   case SdJfast:
      return calcSJ2(energy, redshift);
   case Primack99:
      return calcPrimack(energy, redshift);
   case Primack04:
      return calcPrimack04(energy, redshift);
   }
   return 0;
}

} // namespace IRB
