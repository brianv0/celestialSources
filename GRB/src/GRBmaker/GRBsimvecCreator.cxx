// File: GRBsimvecCreator.cxx

#include <cmath>
#include <algorithm>

#include "GRBsimvecCreator.h"


// static declaration
GRBsimvecCreator	*GRBsimvecCreator::s_instance = 0;




// default constructor
GRBsimvecCreator::GRBsimvecCreator() 
{
    load_dur_loEdge();
    load_dur_long();
    load_dur_short();
    
    load_flux_n();
    load_flux_m();
    load_flux_p();
    load_flux_q();
    
    load_pl_loEdge();
    load_pl_histplaw();
    load_pl_histalpha();
    load_pl_histbeta();
}




// destructor
GRBsimvecCreator::~GRBsimvecCreator ()
{}




// instance()
// Returns the singleton instance of the class
GRBsimvecCreator	*GRBsimvecCreator::instance() 
{ 
    return (s_instance != 0) ? s_instance : (s_instance = new GRBsimvecCreator()); 
}




// kill()
// deletes and releases memory associated with the instance of this class
void GRBsimvecCreator::kill()
{
    delete s_instance;
    s_instance = 0;
}




// load_dur_loEdge(), load_dur_longs(), load_dur_shorts()
// The following three methods load the vectors needed for the calculation of burst durations

void GRBsimvecCreator::load_dur_loEdge()
{
    const double minDur = 0.01;
    const double maxDur = 1000.0;
    
    const long ndecs = long(log10(maxDur)-log10(minDur));
    
    const long binsPerDec = 10;
    
    const long nbins = binsPerDec * ndecs + 1;
    
    const double logStep = 1./binsPerDec;
    
    m_dur_loEdge.resize(nbins);
    for (int i=0; i<nbins; ++i)
        m_dur_loEdge[i] = minDur * pow(10., (logStep*i));
}


void GRBsimvecCreator::load_dur_long()
{
    int dur_long[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,4,0,0,10,20,11,11,
        28,35,44,61,60,83,120,117,98,82,70,53,50,32,22,20,5,9,5,4,0,0,0};
    
    
    long sz = sizeof(dur_long)/sizeof(int);
    
    m_dur_long.resize(sz);
    std::copy(dur_long, dur_long+sz, m_dur_long.begin());
}


void GRBsimvecCreator::load_dur_short()
{
    int dur_short[] = {0,0,0,0,0,0,0,0,7,7,7,27,26,27,34,30,20,27,44,47,57,50,
        24,10,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    
    long sz = sizeof(dur_short)/sizeof(int);
    
    m_dur_short.resize(sz);
    std::copy(dur_short, dur_short+sz, m_dur_short.begin());
}




// load_flux_n(), load_flux_m(), load_flux_p(), load_flux_q()
// The following four methods load the vectors needed for the calculation of peak fluxes

void GRBsimvecCreator::load_flux_n()
{
    long flux_n[] = {8,9,10,12,15,20,25,30,35,40,50,60,70,80,90,100,120,150,200,
        250,300,400,500,600,700,800,900,1000,1100,1200,1262};
    
    long sz = sizeof(flux_n)/sizeof(long);
    
    m_flux_n.resize(sz);
    std::copy(flux_n, flux_n+sz, m_flux_n.begin());
}


void GRBsimvecCreator::load_flux_m()
{
    long flux_m[] = {8,9,10,12,15,20,25,30,35,40,50,60,70,80,90,100,120,150,200,
        250,300,400,430,462};
    
    long sz = sizeof(flux_m)/sizeof(long);
    
    m_flux_m.resize(sz);
    std::copy(flux_m, flux_m+sz, m_flux_m.begin());
}


void GRBsimvecCreator::load_flux_p()
{
    double flux_p[] = {49.279,44.093,42.020,38.082,33.853,28.819,25.417,21.994,
        19.378,18.148,15.910,13.671,12.326,10.387,9.7090,8.5210,7.1180,5.7090,
        4.2970,3.3660,2.7400,1.9810,1.5030,1.2500,1.0350,0.8580,0.7130,0.5950,
        0.5000,0.3890,0.2400};
    
    long sz = sizeof(flux_p)/sizeof(double);
    
    m_flux_p.resize(sz);
    std::copy(flux_p, flux_p+sz, m_flux_p.begin());
}


void GRBsimvecCreator::load_flux_q()
{
    double flux_q[] = {27.912,27.792,24.635,21.766,17.314,15.150,13.985,12.510,
        11.102,9.8450,7.9960,7.3180,6.4150,5.3440,4.7950,4.1910,3.4370,2.9980,
        2.3350,1.9080,1.5720,1.0010,0.8170,0.3710};
    
    long sz = sizeof(flux_q)/sizeof(double);
    
    m_flux_q.resize(sz);
    std::copy(flux_q, flux_q+sz, m_flux_q.begin());
}




// load_pl_loEdge, load_pl_histplaw, load_pl_alpha, load_pl_beta
// The following two methods load the vectors needed in the calculations of power law indices.

void GRBsimvecCreator::load_pl_loEdge()
{
    double value = 7.2;
    
    for (int i=0; i<58; ++i)
        m_pl_loEdge.push_back(value-i*.1);
}


void GRBsimvecCreator::load_pl_histplaw()
{
    int pl_histplaw[] = {1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,1,1,0,0,0,1,0,0,0,0,
        2,0,1,0,0,1,0,3,2,0,0,1,1,0,5,1,4,9,5,8,8,6,6,10,8,3,2,2,1,0};
    
    long sz = sizeof(pl_histplaw)/sizeof(int);
    
    m_pl_histplaw.resize(sz);
    std::copy(pl_histplaw, pl_histplaw+sz, m_pl_histplaw.begin());
}


void GRBsimvecCreator::load_pl_histalpha()
{
    int pl_histalpha[] = {  1,   2,   3,  11,  11,  27,  33, 100, 130, 143,
        227, 270, 364, 390, 426, 508, 435, 398, 365, 271,
        210, 156,  76,  57,  33,  20,  12,  14,  12,   5};
    
    long sz = sizeof(pl_histalpha)/sizeof(int);
    
    m_pl_histalpha.resize(sz);
    std::copy(pl_histalpha, pl_histalpha+sz, m_pl_histalpha.begin());
}


void GRBsimvecCreator::load_pl_histbeta()
{
    int pl_histbeta[] =  { 18,  18,  18,  18,  18,  18,  18,  18,  18,  18,
        18,  18,  18,  18,  18,  18,  18,  18,  18,  18,
        19,  19,  19,  19,  19,  19,  19,  19,  19,  19,
        19,  19,  20,  20,  35,  37,  27,  40,  75, 122,
        117, 160, 167, 215, 263, 310, 333, 365, 450, 487,
        430, 400, 318, 217, 148, 100,  51,  22,   7,   2};
    
    long sz = sizeof(pl_histbeta)/sizeof(int);
    
    m_pl_histbeta.resize(sz);
    std::copy(pl_histbeta, pl_histbeta+sz, m_pl_histbeta.begin());
}


