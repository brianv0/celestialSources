/*!
* \class GRBmaker
*
* \brief This class provides interface for the creation of the bursts of supported types.
* \author Jay Norris        jnorris@lheapop.gsfc.nasa.gov
* \author Sandhia Bansal    sandhiab@lheapop.gsfc.nasa.gov
*
* This class implements methods to generate burst for specific burst parameters as well as to
* generate "n" bursts where "n" is specified in the GRBobsSpectrum.h file.
*
*/


#ifndef GRB_MAKER_H
#define GRB_MAKER_H

#include <vector>
#include <string>

#include "GRBobsUtilities.h"

class GRBurst;
class HepRandomEngine;



class GRBmaker
{
public:
    /*!
     * \brief Constructor.
     */
    GRBmaker()   {}   
    
    /*!
     * \brief Constructor.
     */
    GRBmaker(const double duration, const int npuls, const double flux, 
        const double fraction, const double alpha, const double beta,
        const bool flag);
    
    
    // No memory management function required
    // So - need for destructor, copy constructor and assignment operator to be defined
    
    
    /*!
     * \brief Create "n" bursts in the specified directory
     */
    GRBurst *create(const std::vector<std::string> &paramVector);

    /*!
     * \brief Use specified burst parameters to create photon list
     */
    GRBurst *create(const double duration, const int npuls, const double flux, 
        const double fraction, const double alpha, const double beta, 
        const double epeak, const double specnorm, const bool flag);
    
    /*!
     * \brief Return a copy of itself.
     */
    GRBmaker *clone() const;
    
private:
    // private accessor functions
    
    
    // data members
};


#endif // GRB_MAKER_H
