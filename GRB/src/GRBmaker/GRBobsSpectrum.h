// $Id$

// File: GRBobsSpectrum.h
//
#ifndef GRB_OBS_SPECTRUM_H
#define GRB_OBS_SPECTRUM_H

#include <vector>
#include <string>
#include <functional>  // for unary_function

#include "FluxSvc/ISpectrum.h"

#include "CLHEP/Random/RandomEngine.h"


class GRB;


class GRBobsSpectrum : public ISpectrum
{
    
public:
    GRBobsSpectrum(const std::string &filename);
    GRBobsSpectrum(const double duration, const int npuls, const double flux, const double fraction, const double alpha, 
        const double beta, const double epeak, const double specnorm, const bool flag);
    
    // --- Need a destructor, copy constrctor and assignement operator to manage memory ---
    ~GRBobsSpectrum();
    GRBobsSpectrum(const GRBobsSpectrum &right);   // private copy constructor
    GRBobsSpectrum &operator=(const GRBobsSpectrum &right);
    
    // --- Overridden function ---
    virtual std::string title() const   { return m_title; }
    
    const char * particleName() const   { return m_particleName.c_str(); }
    
    // calculate flux for the current cutoff
    //JCT pure virtual method takes time as argument
    virtual double flux(double time=0) const;
    
    float fraction(float energy);
    
    // return solid angle pair (costh, phi) for the given energy
    std::pair<float,float> dir(float energy) const;
    std::pair<double,double> dir(double energy, HepRandomEngine *engine);
    
    double energySrc(HepRandomEngine *engine, double time=0);
    
    
    
    // --- Functions specific to GRBobsSpectrum class ---
    //virtual float operator() (float randomNumber) const;
    //JCT needs const to match pure virtual method
    virtual float operator() (float randomNumber) const;
    
    // this one asks the GPS for position
    int askGPS();
    
    //JCT pure virtual method 	double ISpectrum::interval(double) missing
    virtual double interval(double time);
    
private:
    
    double nextEnergy() const;
    
    double nextTime() const;
    
    void parseParamList(const std::string &input, std::vector<std::string>& output) const;
    
    void swap(GRBobsSpectrum &other) throw();
    
    //ObserverAdapter< GRBobsSpectrum > m_observer; //obsever tag
    
    
    // std::generate assigns the result of invoking gen, a function object that takes no arguments, 
    // to each element in the range [first, last).
    
    
    struct adder : public std::unary_function<double, double>
    {
        adder(double value) : m_value(value) {}
        double operator() (double x) { return x + m_value; }
        
        double m_value;
    };
    
    
    
    
    
    // Data members:
    std::string        m_title;
    std::string        m_particleName;
    
    GRB               *m_grb;
    
};
#endif // GRB_OBS_SPECTRUM_H

