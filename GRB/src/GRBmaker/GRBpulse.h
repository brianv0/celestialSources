/*!
* \class GRBpulse
*
* \brief This class provides interface for pulse data within a burst.
* \author Jay Norris        jnorris@lheapop.gsfc.nasa.gov
* \author Sandhia Bansal    sandhiab@lheapop.gsfc.nasa.gov
*
*/



#ifndef GRB_PULSE_H
#define GRB_PULSE_H

#include <vector>

class HepRandomEngine;



class GRBpulse
{
public:
    
    /*! 
     * \brief Constructor.
     */
    GRBpulse();
    
    /*! 
     * \brief Destructor.
     */
    ~GRBpulse();

    

    // Accessors
    const std::vector<double> &tmax()  const            { return m_tmax; }
    const std::vector< std::vector<double> > &pulse()   { return m_pulse; }
    const std::vector<long>   &nphotpul() const         { return m_nphotpul; }
    double univFWHM() const                             { return m_univFWHM; }
    


    /*!
     * \brief Return pulse data.
     */
    long data(HepRandomEngine *engine, const double ethres, const long nphoton, 
        const int npuls, const double duration);
    
private:
    
    // Accessors
    const std::vector<double> &amplitude() const        { return m_amplitude; }
    const std::vector<double> &tdiff() const            { return m_tdiff; }
    const std::vector<double> &sigma()  const           { return m_sigma; }

    
    /*!
     * \brief Create vectors needed in the calculations of GRBtimes.
     *
     * See the pulse attributes paper (Norris et al. 1996, ApJ, 459, 393) for description of the 
     * phenomenological "bisigma" pulse model employed here.  
     * Pulse rise-to-decay ratios are ~= 0.4 +- 0.1
     * Generate time profile with 1 millisec precision.  for overall time profile, account for the 
     * T90 duration, the (~ max) rise of 1st pulse (~ 2 s), and the (~ max) decay of last pulse (~ 5 s).
     * Truncate each pulse array @ ~ 5% level, keeping registration with totpulse.
     * Pulses are self-similar:  therefore, number of photons per pulse is proportional to amplitude.
     * Scale the pulse width for energy dependence, as E^-0.333.  make the peak shift 1/2 as large as 
     * width adjustment, ~ as observed at BATSE energies.  
     * These dependences are crude estimates of what we shall measure with GLAST.
     */
    long createSigmaTdiff(HepRandomEngine *engine);
    
    /*!
     * \brief Fill a vector with some hard-coded values.
     */
    void fillVector(std::vector<long> &ngtwid) const;
    
    /*!
     * \brief Return a power-law index from one of the three portions of the pulse width distribution.
     */
    double getAlpha(const double value) const;
    
    /*!
     * \brief Return amplitude for the current burst.
     * 
     * Scramble amplitudes of {1st,2nd} l halves of pulses, separately (leaves whole time profile asymmetric).
     */
    void getAmplitude(HepRandomEngine *engine, const int npuls);
    
    /*!
     * \brief Return number of photons for each pulse in the current burst.
     */
    void getNphotpul(const long nphoton, const int npuls);
    
    /*!
     * \brief Returns pulse data.
     */
    void getPulse(const int npuls);
    
    /*!
     * \brief Return a vector of random numbers multiplied by the duration sorted in ascending order.  
     */
    void getTmax(HepRandomEngine *engine, const int npuls, const double duration);
    
    /*!
     * \brief Return the index i to the last element of "in" vector such that in[i] >= some random value.
     */
    long index(HepRandomEngine *engine, const long diff, const long minval, 
        const std::vector<long> &in) const;
    
    /*! 
     * \brief Choose a universal width for the pulses in the given burst.
     */
    void pickWidth(HepRandomEngine *engine, const double ethres, 
        const double duration);
    
    /*!
     * \brief Calculate the universal width for the pulses in the given burst.
     *
     * Choose a universal width for the pulses within a given burst.
     * A given GRB tends to have pulses of comparable widths.  Therefore (see Fig 3a
     * of Norris et al. 1996 "pulse attributes" paper), pick one pulse width from the
     * distribution of fitted widths of "All" pulses, 50-300 keV, in bright, long
     * BATSE GRBs.  Then, since (a) ~ 1/4 of GRBs are short, and (b) short GRBs have
     * pulse widths ~ 1/10-1/20 that of long GRBs -- multiply pulse widths for one
     * quarter of the GRBs by compression factor of 1/10.  Then using Width ~ E^(-0.333)
     * relationship, scale chosen width (at 100 keV) to width at Ethres.
     */
    void universalWidth(HepRandomEngine *engine, const double ethres, 
        const double duration, const long diff, const long minval, 
        const std::vector<long> &in, const std::vector<double> &v);
    
    
    // Data Members
    std::vector<double>                 m_amplitude;
    std::vector<double>                 m_tdiff;
    std::vector<double>                 m_tmax;
    std::vector<double>                 m_sigma;
    std::vector< std::vector<double> >  m_pulse;
    std::vector<long>                   m_nphotpul;
    double                              m_univFWHM;
};

#endif // GRB_PULSE_H
