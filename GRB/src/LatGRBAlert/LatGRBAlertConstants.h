/*! 
 * \class LatGRBAlertConstants
 * \brief Class instantiated to access general parameters and constants.
 *  
 * The namespace lattrigcst contains all the constant needed by the algorithm to determine
 * trigger likelihoods.
 *
 * The namespace lattrigtypes contains typedef definitions used by the program.
 *
 * \author Jay Norris             jnorris@lheapop.gsfc.nasa.gov
 * \author Sandhia Bansal         sandhiab@lheapop.gsfc.nasa.gov
 */


#ifndef LATGRBALERT_CONSTANTS_H
#define LATGRBALERT_CONSTANTS_H 

#ifndef WIN32
#include <math.h>
#endif



namespace lattrigcst
{
    const double  ratebck     = 4.07;
    const double  degInPi     = 180.0 / M_PI;
    const double  twosigdist  = 34./degInPi;
    
    //	const long    nrange      = 8L;
    const long    nrange      = 80L;
    const long    nmove       = nrange/4;
}



namespace lattrigtypes
{
    typedef std::vector<std::vector<double> >::const_iterator DoubleVecConstIter;
    typedef std::vector<double>::const_iterator DoubleConstIter;
    typedef std::vector<double>::iterator DoubleIter;
    typedef std::vector<double>::const_reverse_iterator DoubleConstRevIter;
}

#endif


