// FILE: LatGrb.cxx


#include <algorithm>              // for sort

#include "LatGrb.h"
#include "GRBpulse.h"
#include "GRBobsConstants.h"
#include "CLHEP/Random/RandFlat.h"


const double  LatGrb::s_emax=100.0;

const double  latGeoArea=282743;
const double  denom=7.;




// Constructor
LatGrb::LatGrb(HepRandomEngine *engine, const std::string &prefix, 
               const std::string &dir) 
      : GRBurst()   
{ 
    createGRB(engine, prefix, dir); 
}



// Constructor
LatGrb::LatGrb(HepRandomEngine *engine, const double duration, const int npuls, 
               const double flux, const double fraction, const double alpha, 
               const double beta, const double epeak, const double specnorm, 
               const bool flag)
      : GRBurst(engine, duration, npuls, flux, fraction, alpha, beta, epeak, 
            specnorm, flag)   
{
}



// calcNphoton(HepRandomEngine *engine)
// Returns number of photons in the current burst
long LatGrb::calcNphoton(HepRandomEngine *engine)
{
    calcSpecnorm(engine);
    
    std::cout << "m_specnorm: " << m_specnorm.back() << " dur: " << 
        m_globalData->duration() << std::endl;
    std::cout << "nphoton: " << 
        long(m_specnorm.back() * m_globalData->duration() * latGeoArea/denom) <<
        std::endl;
    std::cout << "alpha: " << m_globalData->alpha() << " beta: " << 
        m_globalData->beta() << " epeak: " << m_globalData->epeak() << std::endl;
    return long(m_specnorm.back() * m_globalData->duration() * 282743/7.0);
}




// calcSpecnorm(HepRandomEngine *engine)
//		LatGrb::Specnorm computes the spectral normalization for this burst:
//		SpecNorm units, integrated Peak Flux:  photons cm^-2 s^-1 (> Ethres).  Researches on normalization:  
//		(1) J.T. Bonnell's fits to bright BATSE bursts; 
//		(2) comparison with EGRET norms for bright bursts - Catelli's, Dingus' and Schneid's works; and definitively 
//		(3) analysis of Preece et al. (ApJS 2000, 126, 19) spectroscopy catalog of bright BATSE bursts 
//			using IDL routine Specanal.pro.
//
//		The cofactors for N_inc, the number of photons incident on the 6-meter diameter circle are:
//
//		(a) {average flux / peak flux} =~ 1/7;
//		(b) scaling by (peak_flux)^1.5, determined from inspection of Preece et al.;
//		(c) duration (seconds);
//		(d) 282743 cm^2 (6-meter dia. illuminated disk);
//		(e) scaling to integral above Ethres (e.g., 0.03 GeV) for case beta = -2;
//		(f) dispersion (dynamic range) to approximately replicate the scatter in peak flux
//			vs. normalization at 1 MeV as estimated from Preece et al. catalog; and
//      (g) a dependence on power-law index as estimated from Preece et al. catalog.
//
//		Thus, the number of photons produced is consistent with normally incident flux on the projected disk 
//		of the GLAST illumination sphere, integrated above Ethres, for chosen peak flux and duration,  
//		with the photon energies distributed as a power-law.
//
//		Some of these cofactors are applied in LatGrb::specnorm; the remaining ones 
//		are applied in LatGrb::nphoton, which follows this module.
void LatGrb::calcSpecnorm(HepRandomEngine *engine)
{
    // Code for NEW LAT/GBM
    double specnorm = 3.0 * pow((grbcst::ethresLAT*1000.), (-m_globalData->beta() + 1.0)) * pow(m_globalData->fraction(), 1.15);
    specnorm *= (pow(10., (- grbcst::logdyn * engine->flat())));
    std::cout << "specnorm: " << specnorm << std::endl;
    
    
    // Code for OLD LAT
    //double specnorm = 3.0 * pow((grbcst::ethresLAT*1000.), -1) * 
    //    pow(m_globalData->fraction(), 1.5);
    //specnorm *= (25.0 / exp(m_globalData->beta()));
    //specnorm *= (pow(10., (- grbcst::logdyn * engine->flat())));
    
    m_specnorm.push_back(specnorm);
}




// makeEnergies(HepRandomEngine *engine)
//		Choose N_inc energies for this burst, from a single power-law spectral distribution
//		with index = m_beta.  Maximum possible energy to generate is emax, minimum is ethres.
void LatGrb::makeEnergies(HepRandomEngine *engine)
{
    if (m_globalData->beta() == 1.0)
    {
        for (long iphot=0; iphot!=m_nphoton; ++iphot)
        {
            m_photonlist[iphot].setEnergy(grbcst::ethresLAT * exp(engine->flat() * 
                log(s_emax/grbcst::ethresLAT)));
        }
    }
    
    else
    {
        double beta = m_globalData->beta();
        double x = 1.0 - exp((1.0 - beta) * log(s_emax/grbcst::ethresLAT));
        
        for (long iphot=0; iphot!=m_nphoton; ++iphot)
        {
            m_photonlist[iphot].setEnergy(grbcst::ethresLAT * 
                exp(log(1.0 - x*engine->flat()) / (1.0 - beta)));
        }
    }
}



// makeGRB(engine)
// Calls modules which compute energies and times for photons in the current burst.
void LatGrb::makeGRB(HepRandomEngine *engine)
{
    m_photonlist.resize(m_nphoton);
    
    makeEnergies(engine);
    
    makeTimes(engine, grbcst::ethresLAT);
    
    std::sort(m_photonlist.begin(), m_photonlist.end(), timeCmp());
}




