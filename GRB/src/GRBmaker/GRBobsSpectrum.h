/*!
 *\class  GRBobsSpectrum
 *
 * \brief Spectrum class for GRB source physical simulation.
 * Class interfacing the framework with the GRB generation.
 *
 * \author Sandhia Bansal        sandhiab@lheapop.gsfc.nasa.gov
 *
 */

#ifndef GRB_OBS_SPECTRUM_H
#define GRB_OBS_SPECTRUM_H

#include <vector>
#include <string>
#include <functional>  // for unary_function

#include "FluxSvc/ISpectrum.h"

#include "CLHEP/Random/RandomEngine.h"

class GRBurst;



class GRBobsSpectrum : public ISpectrum
{    
public:
    /*! \brief Constructor.
     */
    GRBobsSpectrum(const std::string &filename);

    /*! \brief Constructor.
     */
    GRBobsSpectrum(const double duration, const int npuls, const double flux, 
        const double fraction, const double alpha, const double beta, 
        const double epeak, const double specnorm, const bool flag);
    

    /*! \brief Destructor.
     */
    virtual ~GRBobsSpectrum();
    

    /*! \brief Copy Constructor.
     */
    GRBobsSpectrum(const GRBobsSpectrum &right);   // private copy constructor
    

    /*! \brief Assignment Operator.
     */
    GRBobsSpectrum &operator=(const GRBobsSpectrum &right);

    

    // Accessors
    virtual std::string title() const   { return m_title; }    
    const char * particleName() const   { return m_particleName.c_str(); }
    


    /*! \brief Computes the flux, in \b photons/m^2/s, for a given time
     *  \param time    not used.
     */
    virtual double flux(double time=0) const;
    
    /*! \brief Returns fraction.
     */
    float fraction(float energy);
    
    //! \brief Return galactic direction
    std::pair<float,float> dir(float energy) const;

    //! \brief Calls std::pair<float,float> dir(float energy) const to return galactic direction
    std::pair<double,double> dir(double energy, HepRandomEngine *engine);
    
    /*! \brief returns the energy of a sampled photon.
     *
     *  Method called by \c FluxSource::event(). 
     *  It returns the energy of a sampled photon, by calling 
     *  the \c operator() method. 
     *  \param engine  random engine for uniform sampling;
     *  \param time    current time. 
     */
    double energySrc(HepRandomEngine *engine, double time=0);
       
    /*! \brief Draws from the current spectrum the energy of a sampled photon. 
     *  \param randomNumber Uniform random number drawn in the method \c energySrc .  
     */ 
    virtual float operator() (float randomNumber) const;
    
    /*! \brief Returns the time interval
     */
    virtual double interval(double time);
    
private:
    /*! \brief Returns next energy from the photon list.  Returns an invalid number if at the end of the list
     */
    double nextEnergy() const;
    
    /*! \brief Returns next time from the photon list.  Returns an invalid number if at the end of the list
     */
    double nextTime() const;
    
    /*! \brief Parses input string to extract name of directory and input files and fills output with the values
     */
    void parseParamList(const std::string &input, 
        std::vector<std::string>& output) const;
    
    /*! \brief Used by the assignment operator
     */
    void swap(GRBobsSpectrum &other) throw();
    
    // Data members:
    std::string        m_title;
    std::string        m_particleName;
    GRBurst            *m_grb;   //! Instantiates the GRB model
    
};
#endif // GRB_OBS_SPECTRUM_H

