// File: GRBsimvecCreator.h
// Fills up the vectors needed in the calculations of duration, peak fluxes and power law indices for the GRB bursts

#ifndef GRB_SIMVEC_CREATOR_H
#define GRB_SIMVEC_CREATOR_H

#include <vector>



class GRBsimvecCreator
{
public:
    
    // static access/destruction
    static GRBsimvecCreator	*instance();
    static void     kill ();
    
    // accessor methods:
    const std::vector<double> &dur_loEdge() const    { return m_dur_loEdge; }
    const std::vector<int> &dur_long() const         { return m_dur_long; }
    const std::vector<int> &dur_short() const        { return m_dur_short; }
    const std::vector<long> &flux_n() const          { return m_flux_n; }
    const std::vector<long> &flux_m() const          { return m_flux_m; }
    const std::vector<double> &flux_p() const        { return m_flux_p; }
    const std::vector<double> &flux_q() const        { return m_flux_q; }
    const std::vector<double> &pl_loEdge() const     { return m_pl_loEdge; }
    const std::vector<int> &pl_histplaw() const      { return m_pl_histplaw; }
    const std::vector<int> &pl_histalpha() const     { return m_pl_histalpha; }
    const std::vector<int> &pl_histbeta() const      { return m_pl_histbeta; }
    
protected:
    // singleton - protect ctor/dtor
    GRBsimvecCreator();
    virtual ~GRBsimvecCreator();
    
private:
    static GRBsimvecCreator* s_instance;
    
    std::vector<double>  m_dur_loEdge;
    std::vector<int>     m_dur_long;
    std::vector<int>     m_dur_short;
    
    std::vector<long>    m_flux_n;
    std::vector<long>    m_flux_m;
    std::vector<double>  m_flux_p;
    std::vector<double>  m_flux_q;
    
    std::vector<double>  m_pl_loEdge;
    std::vector<int>     m_pl_histplaw;
    std::vector<int>     m_pl_histalpha;
    std::vector<int>     m_pl_histbeta;
    
    void load_dur_loEdge();
    void load_dur_long();
    void load_dur_short();
    
    void load_flux_n();
    void load_flux_m();
    void load_flux_p();
    void load_flux_q();
    
    void load_pl_loEdge();
    void load_pl_histplaw();
    void load_pl_histalpha();
    void load_pl_histbeta();
};


#endif // GRB_SIMVEC_CREATOR_H
