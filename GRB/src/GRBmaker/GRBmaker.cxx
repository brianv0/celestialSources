// FILE: GRBmaker.cxx

#include "GRBmaker.h"
#include "GRBurst.h"
#include "LatGrb.h"
#include "GbmGrb.h"
#include "CLHEP/Random/RandFlat.h"

using namespace grbobstypes;




// create(const std::vector<std::string> &paramVector)
// Creates LAT and GBM GRBs
GRBurst *GRBmaker::create(const std::vector<std::string> &paramVector)
{
    if (paramVector.size() == 2)
        return new GRBurst(paramVector);
    
    else if (paramVector.size() <= 1)
    {
        // Get an engine and set its seed
        HepRandomEngine *engine = HepRandom::getTheEngine();
        HepRandom::setTheSeed(grbcst::seed);
        
        LatGrb *latGrb = new LatGrb(engine, "LAT", paramVector[0]);
        //GbmGrb *gbmGrb = new GbmGrb(engine, "GBM", latGrb->specnorm(), paramVector[0]);
        
        //delete gbmGrb;
        return latGrb;
    }
    
    else
    {
        std::cout << "invalid number of parameters" << std::endl;
        return 0;
    }
}



// create(const double duration, const int npuls, const double flux, const double fraction, 
//		  const double alpha, const double beta, const double epeak, const double specnorm, const bool flag)
//
// Creates LAT and GBM GRBs
GRBurst *GRBmaker::create(const double duration, const int npuls, const double flux, 
                          const double fraction, const double alpha, 
                          const double beta, const double epeak, 
                          const double specnorm, const bool flag)
{
    // Get an engine and set its seed
    HepRandomEngine *engine = HepRandom::getTheEngine();
    HepRandom::setTheSeed(grbcst::seed);
    
    LatGrb *latGrb = new LatGrb(engine, duration, npuls, flux, fraction, alpha, 
        beta, epeak, specnorm, flag);
    //GbmGrb *gbmGrb = new GbmGrb(engine, duration, npuls, flux, fraction, alpha, 
    //    beta, epeak, specnorm, flag);
    
    //delete gbmGrb;
    return latGrb;
}



// clone()
// Makes a copy of itself and returns it to the caller.
GRBmaker *GRBmaker::clone() const
{
    return new GRBmaker(*this);
}



