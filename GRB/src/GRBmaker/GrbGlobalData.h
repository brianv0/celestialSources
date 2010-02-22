/*!
* \class GrbGlobalData
*
* \brief This class provides interface for the burst specific data.
* \author Jay Norris        jnorris@lheapop.gsfc.nasa.gov
* \author Sandhia Bansal    sandhiab@lheapop.gsfc.nasa.gov
*
* This class follows the SINGLETON pattern.
*/



#ifndef GRB_GLOBAL_DATA_H
#define GRB_GLOBAL_DATA_H

#include <vector>
#include <string>
#include <functional>  // for unary_function

class HepRandomEngine;



class GrbGlobalData
{
public:
    /*!
     * \brief Static access.
     */
    static GrbGlobalData	*instance(HepRandomEngine *engine);

    /*!
     * \brief Static destruction.
     */
    static void            kill ();
    
    
    
    // Accessors
    const std::vector<double> duration() const      { return m_duration; }
    const std::vector<double> flux() const          { return m_flux; }
    const std::vector<double> fraction() const      { return m_fraction; };
    const std::vector<double> alpha() const		    { return m_alpha; }
    const std::vector<double> beta() const			{ return m_beta; }
    const std::vector<double> epeak() const		    { return m_epeak; }
    
    
protected:
    /*!
     * \brief Singleton - protected ctor/dtor
     */
    GrbGlobalData(HepRandomEngine *engine);
    virtual ~GrbGlobalData();
    
    
private:
    /*!
     * \brief instance of the class
     */
    static GrbGlobalData* s_instance;
    
    
    /*!
     * \brief Use interpolation to compute the flux for the burst.
     */
    double computeFlux(HepRandomEngine *engine, const long diff, 
        const long minval, const std::vector<long> &in,
        const std::vector<double> &v) const;
    
    /*!
     * \brief Use interpolation to return the computed value to the caller (duration/power law index).
     */
    double evaluate(HepRandomEngine *engine, const long diff, const long minval,
        const std::vector<long> &in, const std::vector<double> &v) const;
    
    /*!
     * \brief Return samples from global distributions for GRBs: duraton, peak fluxes, and power law indices.
     */
    void grbGlobal(HepRandomEngine *engine);
    
    /*! 
     * \brief Return the index i to the last element of "in" vector such that in[i] < some random value.
     */
    long index(HepRandomEngine *engine, const long diff, const long minval, 
        const std::vector<long> &in) const;
    
    /*!
     * \brief Return burst durations.
     *
     * Chooses durations from the BATSE bimodal duration distribution, where the measurement process is described by 
     * Bonnell et al. (1997, ApJ, 490, 79).  The parent sample is same as for peak fluxes: from GRB 910421 (trig #105)
     * to GRB 990123 (trig #7343).  This partial sample (1262) includes bursts where backgrounds could be fitted, and
     * peak fluxes subsequently measured.  The sample spans 7.75 years.
     */
    void getDuration(HepRandomEngine *engine, long &nlong);
    
    /*!
     * \brief Return fluxes for the burst.
     *
     * Choose peak fluxes from the BATSE log N - log P;  see Bonnell et al. 1997, ApJ, 490, 79, 
     * which duplicates the procedure specified by Pendleton*.  The measurement procedure is applied 
     * uniformly for that part of the BATSE sample from GRB 910421 (trig #105) to GRB 990123 (trig #7343).  
     *
     * This partial sample (1262) includes bursts where backgrounds could be fitted, and peak fluxes subsequently 
     * measured.  It spans 7.75 years.  Therefore, in order to draw from a PF distribution representing 1 year, 
     * we truncate at the eighth brightest burst in 7.75 ~ 8 years.  The peak flux measure in Bonnell et al. is 
     * for 256-ms accumulations.
     *
     * *Pendleton used a different PF estimation technique for the initial BATSE Catalog.
     */
    void getFlux(HepRandomEngine *engine, long nlong);
    
    /*!
     * \brief Choose spectral power-law indices from the BATSE power-law distribution, as measured by Preece et al. (1999).
     */
    void getPowerLawIndex(HepRandomEngine *engine);
    
    /*!
     * \brief Compute power law indices.
     */
    //void powerLawIndex(HepRandomEngine *engine, const std::vector<int> &histpl, 
    //    const std::vector<double> &loEdges,
    //    const double factor, std::vector<double> &vect);
    void powerLawIndex(HepRandomEngine *engine, const std::vector<int> &histpl, const double factor, std::vector<double> &vect);
    
    
    
    // data members
    std::vector<double>   m_duration;
    std::vector<double>   m_flux;
    std::vector<double>   m_fraction;
    std::vector<double>   m_alpha;
    std::vector<double>   m_beta;
    std::vector<double>   m_epeak;
    
    
    /*!
     * \brief Used to fill a vector with the product of its values with m_value.
     */
    struct calcEpeak : public std::unary_function<double, double>
    {
        calcEpeak(HepRandomEngine *engine) : m_engine(engine) {}
        double operator() (double x); 
        
        HepRandomEngine *m_engine;
    };
};



#endif // GRB_GLOBAL_DATA_H
