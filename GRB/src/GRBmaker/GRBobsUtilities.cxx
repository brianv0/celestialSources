// FILE: GRBobsUtilities.cxx

#include <fstream>
#include <iomanip>
#include <algorithm>  // for transform
#include <numeric>  // for accumulate

#include "GRBobsUtilities.h"
#include "GRBsimvecCreator.h"
#include "CLHEP/Random/Random.h"

using namespace grbobstypes;



// result(engine, lo, hi, p)
// returns interpolated values from lo and hi
double GRBobsUtilities::result(HepRandomEngine *engine, const double lo, 
                               const double hi, const double p)
{
    return lo * pow((1.0 - engine->flat() * (1.0 - pow((hi/lo), -p))), (-1/p));
}




// multiplier::operator()
// returns product of x with multiplier::m_value
double GRBobsUtilities::multiplier::operator () (double x)
{
    return x * m_value;
}




// randGen::operator()
// returns a random number between [0,1)
// used to fill a vector with random numbers
double GRBobsUtilities::randGen::operator () (double x)
{
    return m_engine->flat();
}



// sortVector(index, in, sorted, out)
// sorts "out" array in the order of "in"
void GRBobsUtilities::sortVector(const long index, const std::vector<double> &in, 
                                 const std::vector<double> &sorted, 
                                 std::vector<double> &out) 
{
    DoubleSize sz = in.size();
    
    std::vector<double> temp(out);
    
    for (DoubleSize i=0; i<sz; ++i)
    {
        DoubleConstIter it = std::find(sorted.begin(), sorted.end(), in[i]);
        
        DoubleIter it_out = out.begin();
        std::advance(it_out, std::distance(sorted.begin(), it));
        
        *(it_out+index) = temp[i+index];
    }
}