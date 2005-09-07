/**
 * @file EblAtten.h
 * @brief Compute optical depth to extragalactic background light using
 * Hays/McEnery code.
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef eblAtten_EblAtten_h
#define eblAtten_EblAtten_h

#include <map>
#include <stdexcept>

namespace IRB {

/**
 * @class EblAtten
 * @brief Function object wrapper to Hays/McEnery code (in IRB_routines.cxx)
 * that calculates EBL optical depth as a function of energy and redshift
 * for four different models.
 * @author J. Chiang
 *
 * $Header$
 */

enum EblModel {Kneiske, Primack05, Kneiske_HighUV, Salamon_Stecker};

class EblAtten {

public:

   EblAtten(EblModel model);

   /// @return Optical depth to photon-photon absorption.
   /// @param Photon energy (GeV)
   /// @param Source redshift
   float operator()(float energy, float redshift) const;

private:

   EblModel m_model;

   static std::map<EblModel, std::string> s_model_Ids;

};

} // namespace IRB

#endif // eblAtten_EblAtten_h
