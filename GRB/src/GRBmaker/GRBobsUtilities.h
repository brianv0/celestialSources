/*!
* \class GRBobsUtilities
*
* \brief This class provides some utilities for the this program.
* \author Jay Norris        jnorris@lheapop.gsfc.nasa.gov
* \author Sandhia Bansal    sandhiab@lheapop.gsfc.nasa.gov
*
*/


#ifndef GRB_OBS_UTILITIES_H
#define GRB_OBS_UTILITIES_H

#include <vector>
#include <deque>
#include <string>
#include <functional>  // for unary_function and plus

#include "GRBobsConstants.h"

using namespace grbobstypes;

class HepRandomEngine;



class GRBobsUtilities
{
public:
    /*!
     * \breif Return interpolated value from hi and lo.
     */
    static double result(HepRandomEngine *engine, const double lo, 
        const double hi, const double p);
    
    /*!
     * \brief Use ordering of "sorted" to sort "out" vector.
     */
    static void sortVector(const long index, const std::vector<double> &in, 
        const std::vector<double> &sorted, std::vector<double> &out);
     
    /*!
     * \brief Used to fill a vector with the product of its values with m_value.
     */
    struct multiplier : public std::unary_function<double, double>
    {
        multiplier(double value) : m_value(value) {}
        double operator() (double x);
        
        double m_value;
    };
    
    
    /*!
     * =brief Used to fill a vector with random numbers in the range [0,1).
     */
    struct randGen : public std::unary_function<double, double>
    {
        randGen(HepRandomEngine *engine) : m_engine(engine) {}
        double operator() (double x); 
        
        HepRandomEngine *m_engine;
    };
    
    
    /*!
     * \breif Return the ordering of vector "in".
     */
    template<class T>
        static void getSorter(const std::vector<T> &in, 
        std::vector<long> &sorter)
    {
        std::vector<T> sorted(in);
        std::sort(sorted.begin(), sorted.end());
        
        LongSize sz = in.size();
        
        std::deque<bool> qused(sz, 0);
        
        for (LongSize i=0; i<sz; ++i)
        {
            for (LongSize j=0; j<sz; ++j)
            {
                if ((!qused[j]) && (sorted[i] == in[j]))
                {
                    sorter[i] = j;
                    qused[j] = 1;
                    break;
                }
            }
        }
    }
    
    
     /*! 
      * \brief Create cumulative sum of the input vector "in".
      */
    template<class S, class T>
        static void cumulativeSum(const std::vector<S> &in, std::vector<T> &out)
    {
        out.clear();
        
        out.reserve(in.size());
        
        typename std::vector<S>::const_iterator it = in.begin();
        T value = *it;
        out.push_back(value);
        
        ++it;
        
        while (it != in.end())
        {
            value += *it;
            out.push_back(value);
            ++it;
        }
    }

};

#endif
