/*!
* \class GbmGrb
*
* \brief This class is intended to be instantiataed by the GRBmaker::create method.
* \author Jay Norris        jnorris@lheapop.gsfc.nasa.gov
* \author Sandhia Bansal    sandhiab@lheapop.gsfc.nasa.gov
*
*
* This class provides methods to create GBM photon lists in (optionally) specified
* directory.  Alternatively, if the file was generated by a previous run, it can read
* that file and provide data for the simulation.
*/


#ifndef GBM_GRB_H
#define GBM_GRB_H

#include <vector>
#include <string>

#include "GRBurst.h"

class HepRandomEngine;



class GbmGrb  : public GRBurst
{
public:
     /*!
      * \brief Constructor.
      */
    GbmGrb(HepRandomEngine *engine, const std::string &prefix, 
        const std::vector<double> &specnorm, const std::string &dir=0);
    
     /*!
      * \brief Constructor.
      */
    GbmGrb(HepRandomEngine *engine, const double duration, const int npuls, 
        const double flux, const double fraction, const double alpha, 
        const double beta, const double epeak, const double specnorm, 
        const bool flag);
    
     /*!
      * \brief Destructor.
      */
    virtual ~GbmGrb()   {};
    
    

    // Accessors
    static const double emax()   { return s_emax; }
  
    

    
    /*!
    * \brief Calculate number of photons generated by the current burst.
    */
    virtual long calcNphoton(HepRandomEngine *engine, long isim =0);
    
private:
    /*!
    * Choose energies for this burst from a single power-law spectral distribution with 
    * index=lower power law index. 
    */
    void makeEnergies(HepRandomEngine *engine);
    
    /*!
    * Generates n GRBs where n is specified in GRBobsConstants.h file.
    */
    virtual void makeGRB(HepRandomEngine *engine, bool first);
 
    
    
    // Data Members
    static const double  s_emax;
    double               m_cjoin;
    double               m_norm;
};

#endif // GBM_GRB_H
