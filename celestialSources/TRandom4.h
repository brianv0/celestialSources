/**
 * @file TRandom4.h
 * @brief Replace ROOT random number generator with CLHEP version (minimal 
 * version)
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef celestialSources_TRandom4_h
#define celestialSources_TRandom4_h

#include "TRandom3.h"

/**
 * @class TRandom4
 * @brief Replace ROOT uniform random number generator with CLHEP
 * RandFlat.  Do not reimplement sampling of other standard
 * distributions since they appear to rely on the four virtual methods
 * here for generating the needed uniform deviates.
 */

class TRandom4 : public TRandom3 {

public: 
   TRandom4(UInt_t seed=0);
   virtual ~TRandom4();
   virtual Double_t Rndm(Int_t i=0);
   virtual void RndmArray(Int_t n, Float_t *array);
   virtual void RndmArray(Int_t n, Double_t *array);
   virtual void SetSeed(UInt_t seed=0);
};

R__EXTERN TRandom * gRandom;

#endif // celestialSources_TRandom4_h
