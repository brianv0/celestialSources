// FILE: GbmGrb.cxx


#include <algorithm>              // for sort and trasform
#include <iostream>

#include "GbmGrb.h"
#include "GRBpulse.h"
#include "GRBobsConstants.h"


// static constant initialization
const double  GbmGrb::s_emax=10000.0;

// local constant initiaization
const double  gbmGeoArea=126.68;
const double  denom=7.;




// Constructor
GbmGrb::GbmGrb(HepRandomEngine *engine, const std::string &prefix, 
               const std::vector<double> &specnorm, const std::string &dir)
      : GRBurst()
{
    m_specnorm = specnorm;
    createGRB(engine, prefix, dir);
}



// Constructor
GbmGrb::GbmGrb(HepRandomEngine *engine, const double duration, const int npuls, 
               const double flux, const double fraction, const double alpha, 
               const double beta, const double epeak, const double specnorm, 
               const bool flag)
      : GRBurst(engine, duration, npuls, flux, fraction, alpha, beta, epeak, 
      specnorm, flag)
{
}



// calcNphoton(HepRandomEngine *engine)
// Returns number of photons in the current burst.
long GbmGrb::calcNphoton(HepRandomEngine *engine)
{
    static long isim=0;
    
    
    //std::cout << "m_specnorm: " << m_specnorm[isim] << " dur: " << 
    //    m_globalData->duration() << std::endl;
    
    double duration = m_globalData->duration();
    double alpha    = m_globalData->alpha();
    
    double beta     = m_globalData->beta();
    beta = ((beta == 1.0) ? beta+0.001 : beta);
    
    double epeak    = m_globalData->epeak();
    double gamma    = 1.0 - alpha;
    double delta    = 1.0 - beta;
    
    
    
    //std::cout << "duration: " << duration << " alpha: " << 
    //    m_globalData->alpha() << " beta: " << m_globalData->beta() <<
    //    " epeak: " << m_globalData->epeak() << std::endl;
    
    double m_cjoin = pow(epeak, (alpha-beta));
    double m_norm  = m_cjoin * (pow(epeak, gamma) - pow(grbcst::ethresGBM, 
        gamma))/gamma + (pow(s_emax, delta) - pow(epeak, delta))/delta;
    double norm2 = -pow((grbcst::ethresLAT*1.e6), delta) / delta;
    double factBrkplaw = m_norm / norm2;
    
    //std::cout << "m_cjoin: " << m_cjoin << std::endl;
    //std::cout << "gamma: " << gamma << std::endl;
    //std::cout << "delta: " << delta << std::endl;
    //std::cout << "m_norm: " << m_norm << std::endl;
    //std::cout << "norm2: " << norm2 << std::endl;
    //std::cout << "fact: " << factBrkplaw << std::endl;
    //std::cout << "specnorm: " << m_specnorm[isim] << std::endl;
    //std::cout << "nphoton: " << long(factBrkplaw * m_specnorm[isim] * 
    //    duration * gbmGeoArea/denom) << std::endl;
    //exit(0);
    
    return long(factBrkplaw * m_specnorm[isim++] * duration * gbmGeoArea/denom);
}




// makeEnergies(HepRandomEngine *engine)
//		Choose N_inc energies for this burst, from a single power-law spectral distribution
//		with index = m_beta.  Maximum possible energy to generate is emax, minimum is ethres.
void GbmGrb::makeEnergies(HepRandomEngine *engine)
{
    double gamma = 1.0 - m_globalData->alpha();
    double epeak = m_globalData->epeak();
    
    
    std::vector<double> rands(grbcst::nbsim);
    std::transform(rands.begin(), rands.end(), rands.begin(), 
        GRBobsUtilities::randGen(engine));
    
    std::vector<long> dumlo(grbcst::nbsim);
    std::vector<long> dumhi(grbcst::nbsim);
    
    // Fraction of photons below the epeak break
    double fracbreak = m_cjoin * (pow(epeak, gamma) - pow(grbcst::ethresGBM, gamma)) / 
                                                                    (gamma * m_norm);
    
    long nlo=0, nhi=0;
    long i=0;
    for (i=0; i<grbcst::nbsim; ++i)
    {
        if (rands[i] < fracbreak)
            dumlo[nlo++] = i;
        else
            dumhi[nhi++] = i;
    }
    
    dumlo.resize(nlo);
    dumhi.resize(nhi);
    
    double temp;
    long iphot;
    for (iphot=0; iphot<nlo; ++iphot)
    {
        temp = rands[dumlo[iphot]];
        m_photonlist[iphot].setEnergy(exp(log(m_norm*temp*gamma/m_cjoin + 
            pow(grbcst::ethresGBM, gamma)) / gamma));
    }
    
    for (iphot=0; iphot<nhi; ++iphot)
    {
        temp = rands[dumhi[iphot]];
        m_photonlist[iphot].setEnergy(exp(log(m_norm*temp - 
            m_cjoin*(pow(epeak, gamma) - pow(grbcst::ethresGBM, gamma))/gamma)));
    }
}



// makeGRB(engine)
// Calls modules which compute energies and times for photons in the current burst.
void GbmGrb::makeGRB(HepRandomEngine *engine)
{
    m_photonlist.resize(m_nphoton);
    
    makeEnergies(engine);
    
    makeTimes(engine, grbcst::ethresGBM);
    
    std::sort(m_photonlist.begin(), m_photonlist.end(), timeCmp());
}




