// File: GRBpulse.h
// This class creates pulse data for the current burst.

#ifndef GRB_PULSE_H
#define GRB_PULSE_H

#include <vector>


class HepRandomEngine;

class GRBpulse
{
public:
    
    // Constructors/Destructor
    GRBpulse();
    ~GRBpulse();
    
    // Accessor functions
    const std::vector<double> &tmax()  const            { return m_tmax; }
    const std::vector< std::vector<double> > &pulse()   { return m_pulse; }
    const std::vector<long>   &nphotpul() const         { return m_nphotpul; }
    double univFWHM() const                             { return m_univFWHM; }
    
    // Returns pulse data
    long data(HepRandomEngine *engine, const double ethres, const long nphoton, const int npuls, const double duration);
    
    
    
private:
    
    // Private accessor functions
    const std::vector<double> &amplitude() const        { return m_amplitude; }
    const std::vector<double> &tdiff() const            { return m_tdiff; }
    const std::vector<double> &sigma()  const           { return m_sigma; }
    
    // Creates vectors needed in the calculations of GRBtimes
    long createSigmaTdiff(HepRandomEngine *engine);
    
    // Hard codes a vector
    void fillVector(std::vector<long> &ngtwid) const;
    
    // Returns a value based on input value
    double getAlpha(const double value) const;
    
    // Returns amplitude for the current burst
    void getAmplitude(HepRandomEngine *engine, const int npuls);
    
    // Returns number of photons for each pulse in the current burst
    void getNphotpul(const long nphoton, const int npuls);
    
    // pulse data
    void getPulse(const int npuls);
    
    // pulse data
    void getTmax(HepRandomEngine *engine, const int npuls, const double duration);
    
    // returns the index i to the last element of "in" vector such that in[i] >= some random value
    long index(HepRandomEngine *engine, const long diff, const long minval, const std::vector<long> &in) const;
    
    // returns a universal width for the pulses in the given burst
    void pickWidth(HepRandomEngine *engine, const double ethres, const double duration);
    
    // calculates the universal width for the pulses in the given burst
    void universalWidth(HepRandomEngine *engine, const double ethres, const double duration, const long diff, const long minval, 
        const std::vector<long> &in, const std::vector<double> &v);
    
    
    // member data
    std::vector<double> m_amplitude;
    std::vector<double> m_tdiff;
    std::vector<double> m_tmax;
    std::vector<double> m_sigma;
    std::vector< std::vector<double> > m_pulse;
    std::vector<long>   m_nphotpul;
    double m_univFWHM;
};

#endif // GRB_PULSE_H
