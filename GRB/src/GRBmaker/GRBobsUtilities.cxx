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
//
// Input:
//		engine						:	pointer to a HepRandomEngine object
//		lo							:	lower value
//		hi							:	upper value
//		p							:	factor
//
// Output:
//		value						:	interpolated value generated from hi and lo
//
// Calls:
//		engine->flat
//
// Caller:
//		GRBmaker::computeFlux
//		GRBpulse::universalWidth

double GRBobsUtilities::result(HepRandomEngine *engine, const double lo, const double hi, const double p)
{
    return lo * pow((1.0 - engine->flat() * (1.0 - pow((hi/lo), -p))), (-1/p));
}




// multiplier::operator()
// returns product of x with multiplier::m_value
//
// Input:
//		x							:	input parameter
//		m_value						:	multiplier::data member
//
// Output:
//		value						:	x * m_value
//
// Calls:
//		--
//
// Caller:
//		GRBmaker::getTimes
//		GRBpulse::getPulse
//		GRBpulse::getTmax

double GRBobsUtilities::multiplier::operator () (double x)
{
    return x * m_value;
}




// randGen::operator()
// returns a random number between [0,1)
// used to fill a vector with random numbers
//
// Input:
//		--
//
// Output:
//		value						:	random number in [0,1)
//
// Calls:
//		engine->flat
//
// Caller:
//		GRBpulse::getAmplitude
//		GRBpulse::getTmax

double GRBobsUtilities::randGen::operator () (double x)
{
    return m_engine->flat();
}



// sortVector(index, in, sorted, out)
// sorts "out" array in the order of "in"
//
// Input:
//		index						:	offset for the indices into the output array 
//		in							:	input array - unsorted
//		sorted						:	input array - sorted
//		out							:	array to be sorted in the order of "in"
//
// Output:
//		out							:	output array sorted in same order as "in"
//
// Calls:
//		--
//
// Caller:
//		GRBmaker::makeTimes
//		GRBpulse::getAmplitude

// This version is quite slow - but ok for small vectors like the one used by GRBmaker.
void GRBobsUtilities::sortVector(const long index, const std::vector<double> &in, const std::vector<double> &sorted, 
                                 std::vector<double> &out) 
{
    DoubleSize sz = in.size();
    
    std::vector<double> temp(out);
    
    for (long i=0; i<sz; ++i)
    {
        DoubleConstIter it = std::find(sorted.begin(), sorted.end(), in[i]);
        
        DoubleIter it_out = out.begin();
        std::advance(it_out, std::distance(sorted.begin(), it));
        
        *(it_out+index) = temp[i+index];
    }
}