/*!
* \class GbmGrb
*
* \brief This class is intended to be instantiataed by the GRBmaker::create 
*        method.
* \author Sandhia Bansal
*
*
* This class implements methods to create GBM photon lists (time,energy).
*
* If the user specified burst parameters (duration, number of pulses, flux, 
* fraction, broken power law indices, and peak energy), it creates the 
* specific burst and writes it to a file in the current working directory.
*
* If the user specified no input or only a directory name, it creates n GRBs 
* where n is specified in the GRBobsConstants.h file.  If a directory name is 
* specified, the GRBs are written to files in that directory; otherwise, the 
* they are written to files in the current working directory.
*
* The user can also specify names of directory and file containing GRB data in 
* which case, it simply reads the data from the file
*/



#ifndef GBM_GRB_H
#define GBM_GRB_H

#include <vector>
#include <string>

#include "GRB.h"


class HepRandomEngine;



class GbmGrb  : public GRB
{
public:
    /*!
     * \brief Default Constructor.
     */
    GbmGrb(HepRandomEngine *engine, const std::string &prefix, const std::vector<double> &specnorm, const std::string &dir=0);
    
    GbmGrb(HepRandomEngine *engine, const std::string &prefix, const double duration, const int npuls, const double flux,
        const double fraction, const double alpha, const double beta, const double epeak, const double specnorm,
        const bool flag);
    
    // No memory management function required
    
    
    // Accessor Methods
    /*!
     * \brief Returns static constant value.
     */
    static const double emax()   { return s_emax; }
    
    
    // Class Methods
    /*!
     * \brief Calculates number of photons in the current GRB.
     */
    virtual long calcNphoton(HepRandomEngine *engine);
    
private:
    // Class Methods
    /*!
     * Computes photon energies for n GRBs where n is specified in GRBobsConstants.h file.
     */
    void makeEnergies(HepRandomEngine *engine);
    
    /*!
     * Makes n GRBs where n is specified in GRBobsConstants.h file.
     */
    virtual void makeGRB(HepRandomEngine *engine);
    
    
    // Data Members
    static const double  s_emax;
    double               m_cjoin;
    double               m_norm;
};

#endif // GBM_GRB_H
