// FILE: GRBurst.cxx


#include <fstream>
#include <algorithm>              // for transform
#include <numeric>                // for accumulate
#ifdef WIN32
#include <sstream>
#else // gcc < 2.9.5.3 dosen't know <sstream>,stringstream
#include <sstream>
#include <strstream>
#endif

#include "GRBurst.h"
#include "GRBpulse.h"
#include "GrbGlobalData.h"
#include "GRBobsConstants.h"
#include "GRBsimvecCreator.h"
#include "CLHEP/Random/RandFlat.h"
#include "facilities/Util.h"


using namespace grbobstypes;



// Constructor  GRBurst()
// Initializes member data
GRBurst::GRBurst()
       : m_grbdir(std::make_pair<float,float>(0.0,0.0)),
         m_univFWHM(0.0),
         m_specnorm(),
         m_nphoton(0),
         m_photonlist(),
         m_globalData(new GlobalData)
{
}



// Constructor  GRBurst(const std::vector<std::string> &paramVector)
// Initializes member data and reads GRBurst from specified file in paramVector
GRBurst::GRBurst(const std::vector<std::string> &paramVector)
       : m_grbdir(std::make_pair<float,float>(0.0,0.0)),
         m_univFWHM(0.0),
         m_specnorm(),
         m_nphoton(0), 
         m_photonlist(),
         m_globalData(new GlobalData)
{
    readGRB(paramVector);
}



GRBurst::GRBurst(HepRandomEngine *engine, const double duration, const int npuls, 
         const double flux, const double fraction, const double alpha, 
         const double beta, const double epeak, const double specnorm, 
         const bool flag)
       : m_grbdir(std::make_pair<float,float>(0.0,0.0)),
         m_univFWHM(0.0),
         m_specnorm(),
         m_nphoton(0),
         m_photonlist(),
         m_globalData(new GlobalData)
{
    createGRB(engine, duration, npuls, flux, fraction, alpha, beta, epeak, 
        specnorm, flag);
}



// Copy Constructor
// Initializes member data and reads GRBurst from specified file in paramVector
GRBurst::GRBurst(const GRBurst &right)
       : m_grbdir(right.m_grbdir),
         m_univFWHM(right.m_univFWHM),  
         m_specnorm(right.m_specnorm),
         m_nphoton(right.m_nphoton),
         m_photonlist(right.m_photonlist)
{
    if (right.m_globalData)
        m_globalData = right.m_globalData->clone();
}



// Destructor
GRBurst::~GRBurst()
{
    delete m_globalData;
    m_globalData = 0;
}



// Helper function for the assignment operator
void GRBurst::swap(GRBurst &other) throw()
{
    std::swap(m_grbdir, other.m_grbdir);
    std::swap(m_univFWHM, other.m_univFWHM);
    std::swap(m_specnorm, other.m_specnorm);
    std::swap(m_nphoton, other.m_nphoton);
    std::swap(m_photonlist, other.m_photonlist);
    std::swap(m_globalData, other.m_globalData);
}



// Assignment Operator
GRBurst &GRBurst::operator=(const GRBurst &right)
{
    GRBurst temp(right);   // does all the work
    swap(temp);   // this can't throw
    return *this;
}


// createGRB(HepRandomEngine *engine, const std::string &prefix, const std::string &dir)
//		GRB Simulation Overview:
//		(1) GRBurst calls GRBglobal, which calls procedures that return Nbsim samples from each of the burst
//			integral distributions { durations, peak fluxes, power-law indices }.
//		(2) Module makeGRB computes the number of photons for this burst.
//		(3) Module maketimes makes BATSE-like time profiles, but with pulse widths extrapolated to
//			GLAST energies, generated in a call to module pickwidth; photons are distributed within a 
//			given pulse according to an energy-dependent formulation.
//
//		Generates "nbsim" number of bursts.  Loops through each burst, creates a photon list (time,energy) 
//      and records the list in an output file.
void GRBurst::createGRB(HepRandomEngine *engine, const std::string &prefix, 
                        const std::string &dir)
{
    // Choose from the distributions for durations, peak fluxes, and spectral power-law indices
    std::vector<double> duration  = GrbGlobalData::instance(engine)->duration();	
    std::vector<double> flux      = GrbGlobalData::instance(engine)->flux();	
    std::vector<double> fraction  = GrbGlobalData::instance(engine)->fraction();	
    std::vector<double> alpha     = GrbGlobalData::instance(engine)->alpha();	
    std::vector<double> beta      = GrbGlobalData::instance(engine)->beta();	
    std::vector<double> epeak	  = GrbGlobalData::instance(engine)->epeak();	
    
    // create base name for the output files to be generated
    std::string base = baseFilename(prefix, dir);
    
    
    // For each burst, create its attributes
    for (long isim=0; isim<grbcst::nbsim; ++isim)
    {
        std::cout << "Processing " << prefix << " Burst: " << isim << std::endl;
        //std::cout << "dur: " << duration[isim] << " alpha: " << 
        //    alpha[isim] << " beta: " << beta[isim] <<
        //    " epeak: " << epeak[isim] << std::endl;
        
        m_globalData->setDuration(duration[isim]);
        m_globalData->setFlux(flux[isim]);
        m_globalData->setFraction(fraction[isim]);
        m_globalData->setAlpha(alpha[isim]);
        m_globalData->setBeta(beta[isim]);
        m_globalData->setEpeak(epeak[isim]);
        m_globalData->setNpuls(0);
        
        //std::cout << "call calcNphoton\n";
        if ((m_nphoton=calcNphoton(engine, isim)) > 5)
        {
            m_globalData->setSpecnorm(m_specnorm[isim]);
            m_grbdir = direction(engine);
            
            makeGRB(engine, (isim==0));
            
            std::string fname = outputFilename(base, isim);
            std::ofstream os(fname.c_str());
            os << *this;
        }
        
        else
        {
            std::cout << "no GRBurst made: Insufficient flux;" << std::endl;
            std::cout << "nphoton: " << m_nphoton << std::endl << std::endl;
        }
    }
    
    //	GrbGlobalData::kill();
}




// createGRB(HepRandomEngine *engine, const std::string &prefix, const double duration, const int npuls, const double flux,
//			   const double fraction, const double alpha, const double beta, const double epeak, const double specnorm,
//			   const bool flag)
// Creates the photon list (time,energy) for the burst specified by the input parameters 
//		and records it in a file if the flag is set.
void GRBurst::createGRB(HepRandomEngine *engine, const double duration, 
                        const int npuls, const double flux, const double fraction, 
                        const double alpha, const double beta, const double epeak, 
                        const double specnorm, const bool flag)
{
    m_globalData->setDuration(duration);
    m_globalData->setFlux(flux);
    m_globalData->setFraction(fraction);
    m_globalData->setAlpha(alpha);
    m_globalData->setBeta(beta);
    m_globalData->setEpeak(epeak);
    m_globalData->setSpecnorm(specnorm);
    m_globalData->setNpuls(npuls);
    
    if ((m_nphoton=calcNphoton(engine)) > 5)
    {
        m_grbdir = direction(engine);
        
        std::cout << "Creating photon list..." << std::endl;
        makeGRB(engine, true);
        
        if (flag)
        {
            std::string out("GRB_c3.lis");
            std::ofstream os(out.c_str());
            os << *this;
        }
    }
    
    else 
    {
        std::cout << "no GRBurst made: Insufficient flux;" << std::endl;
        std::cout << "nphoton: " << m_nphoton << std::endl;
    }
}




// readGRB(const std::vector<std::string> &paramVector)
// Reads the photon list (time,energy) generated by an instantiation of GRBurst().
void GRBurst::readGRB(const std::vector<std::string> &paramVector)
{
    try
    {
        std::string location = paramVector[0];
        //std::cout << "location: " << location << std::endl;
        facilities::Util::expandEnvVar(&location);
        
        std::string fname = location + paramVector[1];
        
        std::cout << "Reading photon list from " << fname << "..." << std::endl;
        std::ifstream is (fname.c_str());
        is >> *this;
        
        //std::cout << "m_flux: " << m_globalData->flux() << std::endl;
        std::ofstream os("GRB_c2.lis");
        os << *this;
    }
    
    catch (...)
    {
        std::cout << "Error while opening input file: " << paramVector[1] << std::endl;
    }
}




// clone()
// Makes a copy of itself and returns it to the caller.
GRBurst *GRBurst::clone() const
{
    return new GRBurst(*this);
}



// Default constructor  baseFilename(const std::string &prefix, const std::string &dir)
// Returns a string containing base name used for the generation of output file name.
std::string GRBurst::baseFilename(const std::string &prefix, const std::string &dir) const
{
    std::string base;
    
    if (dir.size() > 1)    // null terminated - so the size is at least "1"
    {
        std::string location = dir;
        facilities::Util::expandEnvVar(&location);
        base = location + prefix + "_GRB_";
    }
    
    else
        base = prefix + "_GRB_";
    
#ifdef WIN32    
    std::ostringstream ind;
    ind << grbcst::nbsim;
    int i = strlen(ind.str().c_str());
#else
    std::ostringstream ind;
    ind << grbcst::nbsim;
    int i = strlen(ind.str().c_str());
    //    char *ind = new char(80);
    //    sprintf(ind, "%ld", grbcst::nbsim);
    //    int i = strlen(ind);
#endif
    
    for (int j=0; j<i; ++j)
        base += '0';
    
    return base;
}




// Default constructor  outputFilename()
// Uses base string generated by baseFilename to create name of the output file.
std::string GRBurst::outputFilename(const std::string &base, const long isim) const
{
    std::string::size_type baseSize = base.size();
    
    std::string fname = base;
    
#ifdef WIN32
    std::ostringstream ind;
    ind << isim;
    fname.replace(baseSize-strlen(ind.str().c_str()), baseSize-1, ind.str());
#else
    std::ostringstream ind;
    ind << isim;
    fname.replace(baseSize-strlen(ind.str().c_str()), baseSize-1, ind.str());
    //    char *ind = new char(80);
    //    sprintf(ind, "%ld", isim);
    //    fname.replace(baseSize-strlen(ind), baseSize-1, ind);
#endif
    
    fname += ".lis";
    
    return fname;
}




// direction(engine)
// Calculates the direction of the burst.
std::pair<float,float> GRBurst::direction(HepRandomEngine *engine) const
{
    float coszenith = 1.0 - engine->flat() * grbcst::zenNorm;
    float azimuth    = 2. * M_PI * engine->flat();
    
    return std::make_pair<float,float> (-coszenith, azimuth);
}




// getTimes(HepRandomEngine *engine, const double ethres, const long nphots, const long deltbinsleft, const long iphotoff, 
//				   const double tmax, const std::vector<double> &pulse)
// Uses pulse data to generate list of photon times
void GRBurst::getTimes(HepRandomEngine *engine, const double ethres, 
                       const long nphots, const long deltbinsleft, 
                       const long iphotoff, const double tmax, 
                       const std::vector<double> &pulse)
{
    std::vector<double> cumpulse;
    
    GRBobsUtilities::cumulativeSum(pulse, cumpulse);
    double maxval = cumpulse[cumpulse.size()-1];
    std::transform(cumpulse.begin(), cumpulse.end(), cumpulse.begin(), 
        GRBobsUtilities::multiplier(1./maxval));
    
    double photontime, efactor, twidthscale;
    
    for (int iphot=0; iphot<nphots; ++iphot)
    {
        DoubleIter it = cumpulse.end();
        while (it == cumpulse.end())
            it = std::find_if(cumpulse.begin(), cumpulse.end(), 
                std::bind2nd(std::greater_equal<double>(), engine->flat()));
        
        int index = std::distance(cumpulse.begin(), it);
        
        photontime   = (index - deltbinsleft) * grbcst::timres;
        efactor      = pow((m_photonlist[iphot+iphotoff].energy() / ethres), -0.333);
        
        twidthscale  = photontime * efactor;
        
        m_photonlist[iphot+iphotoff].setTime(twidthscale + tmax + 
            0.5 * efactor * deltbinsleft * grbcst::timres);
    }
}




// makeTimes(HepRandomEngine *engine, bool first, double ethres)
//		Makes BATSE-like GRBurst time profiles, placing GLAST photons a la cumulative BATSE intensity, 
//		but in narrower pulses:
//		(1) m_npuls = number of pulses, proportional to BATSE duration.
//		(2) pulse peak amplitude is random (0.0=>1.0); sort amps in descending amp order.
//		(3)	scramble amps of {1st,2nd} halves of pulses, separately (leaves profile assymmetric).
//		(4)	center of pulse time is random within duration.  sort the times in ascending order.
//		(5)	pulse width is drawn from BATSE width distribution for bright bursts (attributes paper), 
//			scaled to GLAST energies using width ~E^-0.333.
//		(6)	make m_npuls pulses with "bisigma" shapes => sum to produce time profile
//		(7)	form cumulative distribution of BATSE-like intensity
//		(8) distribute the m_nphoton photons according to cumulative intensity => GRBtimes
//		(9)	offset the photon times according to
//			(a) energy dependence, width ~E^-0.333 and
//			(b) time of peak, also proportional to E^-0.333.
void GRBurst::makeTimes(HepRandomEngine *engine, bool first, double ethres)
{
    if (m_globalData->npuls() == 0)
        m_globalData->setNpuls(std::max(int(15 * m_globalData->duration()/30 + 0.5), 1));
    
    try
    {
        GRBpulse *grbPulse = new GRBpulse;
        int npuls = m_globalData->npuls();
        
        long deltbinsleft = grbPulse->data(engine, first, ethres, m_nphoton, npuls, 
            m_globalData->duration());
        
        // create times for each pulse
        const std::vector<double> tmax                 = grbPulse->tmax();
        const std::vector< std::vector<double> > pulse = grbPulse->pulse();
        const std::vector<long> nphotpul               = grbPulse->nphotpul();
        
        long iphotoff = 0;
        for (long ipuls=0; ipuls<npuls; ++ipuls)
        {
            if (nphotpul[ipuls] > 0)
            {
                getTimes(engine, ethres, nphotpul[ipuls], deltbinsleft, 
                    iphotoff, tmax[ipuls], pulse[ipuls]);
                iphotoff += nphotpul[ipuls];
            }
        }
        
        m_univFWHM = grbPulse->univFWHM();
        
        delete grbPulse;
    }
    
    catch(...)
    {
        std::cout << "Cannot allocate memory for the GRBurst pulse..." << std::endl;
    }
}




// operator>>(is, grb)
// Reads burst information and photon list (time,energy) from the input stream and updates the object grb with it
std::ifstream &operator>>(std::ifstream &is, GRBurst &grb)
{
    std::string str;
	// >>>>>>>>>>>>>> TEMPORARILY COMMENT OUT READING OF HEADER FOR NOW -- 12/08/03 <<<<<<<<<<<<<<<
    //is >> str >> str >> str >> str >> str >> str;
    
    double duration, flux, beta, specnorm;
    int    npuls;
    
	// >>>>>>>>>>>>>> TEMPORARILY COMMENT OUT READING OF HEADER FOR NOW -- 12/08/03 <<<<<<<<<<<<<<<
    //is >> flux >> duration >> beta >> specnorm >> npuls >> grb.m_univFWHM;
    grb.m_globalData->setFlux(2.5);
    //grb.m_globalData->setFlux(flux);
    //grb.m_globalData->setFlux(flux);
    //grb.m_globalData->setDuration(duration);
    //grb.m_globalData->setBeta(beta);
    //grb.m_globalData->setNpuls(npuls);
    //grb.m_globalData->setSpecnorm(specnorm);

	// >>>>>>>>>>>>>> TEMPORARILY COMMENT OUT READING OF ZENITH/AZIMUTH FOR NOW -- 12/08/03 <<<<<<<<<<<<<<<
    //is >> str >> str >> str;
    //float zenith, azimuth;
    //is >> zenith >> azimuth;
    //grb.m_grbdir = std::make_pair<float,float> (zenith, azimuth);
    
    //is >> str >> str;
    is >> grb.m_nphoton;
    //is >> str >> str >> str >> str;
    
    grb.m_photonlist.resize(grb.m_nphoton);
    
    //double t, e;
	double t, e, ra, dec;
    for (long i=0; i<grb.m_nphoton; ++i)
    {
        //s >> t >> e;

		// >>>>>>>>>>>>>>> TEMPORARILY READ TIME, ENERGY, RA, DEC FROM PHOTON LIST -- 12/08/03 <<<<<<<<<<<<<<<<<<<
        is >> t >> e >> ra >> dec;
        
        grb.m_photonlist[i].setTime(t);
        grb.m_photonlist[i].setEnergy(e);

	std::pair<float,float> radec(ra,dec);
	grb.m_photonlist[i].setDir(radec);
    }
    
	std::cout << "Finished reading file...\n";
    return is;
}




// operator>>(os, grb)
// Writes burst information and photon list (time,energy) contained in the object grb to the output stream
std::ofstream &operator<<(std::ofstream &os, const GRBurst &grb)
{
    std::pair<double,double> grbdir = grb.dir();
    
    GlobalData *globalData = grb.globalData();
    
    os << "Fp, dur, beta, Specnorm, Npulses, UnivFWHM= " << std::setw(12) << 
        std::setiosflags(std::ios::fixed) << 
        globalData->flux() << "   " <<
        globalData->duration() << "   " << globalData->beta() << "   " << 
        globalData->specnorm() << "   " << 
        globalData->npuls() << "  " << grb.univFWHM() << std::endl;
    
    os << "ZenAng, AziAng = " << std::setw(12) << 
        std::setiosflags(std::ios::fixed) << grbdir.first << "   " << 
        grbdir.second << std::endl;
    
    os << "nphotons = " << grb.nphoton() << std::endl;
    
    os << "(Times, Energies) per photon:" << std::endl;
    
    
    long nphoton = grb.nphoton();
    if (nphoton > 0)
    {
         std::vector<TimeEnergy> photonlist = grb.photonlist();
        for (long i=0; i<nphoton; ++i)
        {
            double time = photonlist[i].time();
            os << std::setw(12) << std::setiosflags(std::ios::fixed) << 
                photonlist[i].time() << "   " << 
                photonlist[i].energy() << std::endl;
        }
    }
    
    return os;
}
