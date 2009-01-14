/**
 * @file EblAtten.h
 * @brief Compute optical depth to extragalactic background light using
 * Hays/McEnery code.
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef IRB_EblAtten_h
#define IRB_EblAtten_h

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

enum EblModel {SdJbase, SdJfast, Primack99, Primack04};

class EblAtten {

public:

   EblAtten(EblModel model);

   /// @return Optical depth to photon-photon absorption.
   /// @param Photon energy (MeV)
   /// @param Source redshift
   float operator()(float energy, float redshift) const;

private:

   EblModel m_model;

};

} // namespace IRB

#endif // IRB_EblAtten_h
