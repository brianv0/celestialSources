// FILE: GRBpulse.cxx

#include <algorithm>  // for transform
#include <numeric>  // for accumulate

#include "GRBpulse.h"
#include "GRBobsUtilities.h"
#include "GRBobsConstants.h"
#include "CLHEP/Random/Random.h"

using namespace grbobstypes;



// Default constructor  GRBpulse()
// Initializes the member data.
//
// Input:
//		--
//
// Output:
//		GRBpulse object
//
// Calls:
//		--

GRBpulse::GRBpulse()
:m_amplitude(),
m_tdiff(),
m_tmax(),
m_sigma(),
m_pulse(),
m_nphotpul(),
m_univFWHM(0.0)
{
}




// Default destructor ~GRBpulse()
GRBpulse::~GRBpulse()
{
}




// getAmplitude(HepRandomEngine *engine, const int npuls)
//		Generates pulse amplitudes and sorts them in ascending order.
//		Scramples amplitudes of {1st,2nd} l halves of pulses, separately (leaves whole time profile asymmetric)
//
// Input:
//		engine						:	pointer to a HepRandomEngine object
//		npuls				 		:	number of pulses in this burst
//
// Output:
//		m_amplitude					:	pulse peak amplitude
//
// Calls:
//		GRBobsUtilities::sortVector
//		GRBobsUtilities::randGen
//
// Caller:
//		data

void GRBpulse::getAmplitude(HepRandomEngine *engine, const int npuls)
{
    m_amplitude.resize(npuls);
    std::transform(m_amplitude.begin(), m_amplitude.end(), m_amplitude.begin(), GRBobsUtilities::randGen(engine));
    std::sort(m_amplitude.begin(), m_amplitude.end(), std::greater<double>());
    
    if (npuls > 3)
    {
        long n1st = npuls / 2;
        long n2nd = n1st;
        if (npuls % 2 == 1)
            ++n1st;
        
        std::vector<double> index1st(n1st);
        std::transform(index1st.begin(), index1st.end(), index1st.begin(), GRBobsUtilities::randGen(engine));
        
        std::vector<double> sorted(index1st);
        std::sort(sorted.begin(), sorted.end());
        GRBobsUtilities::sortVector(0, index1st, sorted, m_amplitude);
        
        index1st.resize(n2nd);
        std::transform(index1st.begin(), index1st.end(), index1st.begin(), GRBobsUtilities::randGen(engine));
        
        sorted = index1st;
        std::sort(sorted.begin(), sorted.end());
        GRBobsUtilities::sortVector(n1st, index1st, sorted, m_amplitude);
    }
}




// getAlpha(const double value)
//		Gets a power-law index from one of the three portions of the pulse width distribution
//
// Input:
//		value						:	input value
//
// Output:
//		value						:	output value
//
// Calls:
//		---
//
// Caller:
//		universalWidth

double GRBpulse::getAlpha(const double value) const
{
    if (value >= 0.01 && value <= 0.10)
        return -0.1;
    
    if (value > 0.10 && value <= 0.25)
        return -1.5;
    
    if (value > 0.25 && value <= 0.40)
        return 0.25;
    
    if (value > 0.40)
        return 0.6;
    
    return 0;
}




// fillVector(std::vector<long> &ngtwid)
//		Fills up a vector with the hard coded values for the integral pulse width distribution.
//
// Input:
//		---
//
// Output:
//		ngtwid						:	filled vector
//
// Calls:
//		---
//
// Caller:
//		pickWidth

void GRBpulse::fillVector(std::vector<long> &ngtwid) const
{
    // Keep it in sorted order for index method to work
    // It can be sorted by the caller - however, it is just as easy to keep it as sorted
    
    long temp[] = {430,427,424,421,417,410,388,334,248,178,119,81,46,15,2};
    
    LongSize sz = sizeof(temp)/sizeof(long);
    
    ngtwid.resize(sz);
    std::copy(&temp[0], &temp[sz], ngtwid.begin());
}




// index(HepRandomEngine *engine, const long diff, const long minval, const std::vector<long> &in)
// returns the index i to the last element of "in" vector such that in[i] >= some random value.
//
// Input:
//		engine						:	pointer to a HepRandomEngine object
//		in							:	vector to be scanned for some specified value - SORTED IN ASCENDING ORDER
//		diff						:	difference between maximum and minimum elements of vector "in"
//		minval						:	minimum element of vector "in"
//
// Output:
//		index						:	index such that in[index] >= some random value
//
// Calls:
//		engine->flat
//
// Caller:
//		universalWidth

long GRBpulse::index(HepRandomEngine *engine, const long diff, const long minval, const std::vector<long> &in) const
{
    bool found=0;
    long i;
    static long sz = in.size();
    
    while (!found)
    {
        // use random number to generate a threshold value
        const long value = long(engine->flat() * diff + minval);
        
        // pick last index i into "in" such that in[i] < value
        for (i=in.size()-1; i>=0 && !found; --i)
            if (value <= in[i])
                found = 1;
    }
    
    return i+1;
}




// universalWidth(HepRandomEngine *engine, const double ethres, const double duration, const long diff, 
//				  const long minval, const std::vector<long> &in, const std::vector<double> &v)
//		Chooses a universal width for the pulses within a given burst.
//		A given GRB tends to have pulses of comparable widths.  Therefore (see Fig 3a
//		of Norris et al. 1996 "pulse attributes" paper), pick one pulse width from the
//		distribution of fitted widths of "All" pulses, 50-300 keV, in bright, long
//		BATSE GRBs.  Then, since (a) ~ 1/4 of GRBs are short, and (b) short GRBs have
//		pulse widths ~ 1/10-1/20 that of long GRBs -- multiply pulse widths for one
//		quarter of the GRBs by compression factor of 1/10.  Then using Width ~ E^(-0.333)
//		relationship, scale chosen width (at 100 keV) to width at Ethres.
//
// Input:
//		grbcst::ethres				:	double
//		engine						:	pointer to a HepRandomEngine object
//		duration					:	current burst duration
//		v							:	hi and lo values picked from vector "v" are interpolated to get the return value
//		in							:	vector to aid in picking values from the "v" vector
//		diff						:	difference between maximum and minimum elements of vector "in"
//		minval						:	minimum element of vector "in"
//
// Output:
//		ngtwid						:	filled vector
//
// Calls:
//		index
//		GRBobsUtilities::result
//
// Caller:
//		pickWidth

void GRBpulse::universalWidth(HepRandomEngine *engine, const double ethres, const double duration, const long diff, 
                              const long minval, const std::vector<long> &in, const std::vector<double> &v) 
{
    // find index loIndex such that in[loIndex] < some random number
    long loIndex = index(engine, diff, minval, in);
    
    double wid_lo = v[loIndex];
    double wid_hi = v[loIndex+1];
    
    double p = getAlpha(wid_lo);
    
    m_univFWHM = GRBobsUtilities::result(engine, wid_lo, wid_hi, p);
    
    if (duration < 2.5)
        m_univFWHM /= 10;
    
    m_univFWHM *= pow((ethres/0.0001), -0.333);
}




// pickWidth(HepRandomEngine *engine, const double ethres, const double duration)
// Chooses a universal width for the pulses within a given burst.
//
// Input:
//		grbcst::nbins				:	long
//		grbcst::logfac0				:	double
//		engine						:	pointer to a HepRandomEngine object
//		duration					:	current burst duration
//
// Output:
//		m_univFWHM					:	universal width for the pulses within current burst
//
// Calls:
//		fillVector
//		universalWidth
//
// Caller:
//		data

void GRBpulse::pickWidth(HepRandomEngine *engine, const double ethres, const double duration)
{
    static std::vector<double> logwidth;
    static std::vector<long> ngtwid;
    static long minNgtwid;
    static long diff;
    
    if (ngtwid.empty())
    {
        int sz = grbcst::nbins + 1;
        
        logwidth.reserve(sz);
        
        for (int i=0; i<sz; ++i)
            logwidth.push_back(grbcst::minwid * pow(grbcst::logfac0, i));
        
        fillVector(ngtwid);
        
        LongConstIter it = std::max_element(ngtwid.begin(), ngtwid.end());
        long maxNgtwid = std::max(*it, 0L);
        
        it = std::min_element(ngtwid.begin(), ngtwid.end());
        minNgtwid = std::min(*it, 0L);
        
        diff = maxNgtwid - minNgtwid;
    }
    
    universalWidth(engine, ethres, duration, diff, minNgtwid, ngtwid, logwidth);
}




// createSigmaTdiff(HepRandomEngine *engine)
//		See the pulse attributes paper (Norris et al. 1996, ApJ, 459, 393) for
//		description of the phenomenological "bisigma" pulse model employed here.
//		pulse rise-to-decay ratios are ~= 0.4 +- 0.1
//		generate time profile with 1 millisec precision.  for overall time profile,
//		account for the T90 duration, the (~ max) rise of 1st pulse (~ 2 s), and
//		the (~ max) decay of last pulse (~ 5 s).
//		truncate each pulse array @ ~ 5% level, keeping registration with totpulse.
//		pulses are self-similar:  therefore, number of photons per pulse is proportional to amplitude
//		scale the pulse width for energy dependence, as E^-0.333.  make the peak shift
//		1/2 as large as width adjustment, ~ as observed at BATSE energies.
//		theses dependences are crude estimates of what we shall measure with GLAST.

//
// Creates the two vectors sigma and tdiff needed in the calculation of GRB times.
//
// Input:
//		grbcst::nu					:	long
//		grbcst::frackeep			:	double
//		grbcst::timres				:	double
//		engine						:	pointer to a HepRandomEngine object
//		m_univFWHM					:	universal width for the pulses within the current burst
//
// Output:
//		deltbinsleft				:	value used in the calculations of GRB times
//		m_sigma						:	vector containing values equal to sigrise and sigdecay
//		m_tdiff						:	
//
// Calls:
//		engine->flat
//
// Caller:
//		data

long GRBpulse::createSigmaTdiff(HepRandomEngine *engine)
{
    double riseHWHM  = (0.4 + 0.1 * engine->flat()) * m_univFWHM;
    double decayHWHM = (1.-riseHWHM) * m_univFWHM;
    
    double denom     = pow(log(2.), (1.0/grbcst::nu));
    double sigrise   = riseHWHM / denom;
    double sigdecay  = decayHWHM / denom;
    
    double factor         = pow(-log(grbcst::frackeep), (1.0/grbcst::nu));
    long   deltbinsleft   = long((sigrise * factor) / grbcst::timres + 1);
    long   deltbinsrite   = long((sigdecay * factor) / grbcst::timres + 1);
    long   lenrange       = deltbinsleft + deltbinsrite;
    
    
    m_sigma.resize(lenrange, sigrise);
    m_tdiff.reserve(lenrange);
    
    double prod = deltbinsleft * grbcst::timres;
    double value;
    
    for (long i=0; i!=lenrange; ++i)
    {
        m_tdiff.push_back(value=(i * grbcst::timres - prod));
        if (value >= 0)
            m_sigma[i] = sigdecay;
    }
    
    return deltbinsleft;
}




// getNphotpul(const long nphoton, const int npuls)
// Computes number of photons in each pulse for the given burst.
//
// Input:
//		nphoton						:	number of photons in the current burst
//		npuls						:	number of pulses in the current burst
//		m_amplitude					:	peak flux amplitude
//
// Output:
//		m_nphotpul					:	vector containing number of photons in each pulse
//
// Calls:
//		GRBobsUtilities::getSorter
//
// Caller:
//		data

void GRBpulse::getNphotpul(const long nphoton, const int npuls)
{
    long v;
    long totnpht = 0;
    double totamp = std::accumulate(m_amplitude.begin(), m_amplitude.end(), 0.0);
    
    m_nphotpul.resize(npuls);
    for (int i=0; i<npuls; ++i)
    {
        v = long((m_amplitude[i]/totamp) * nphoton + 0.5);
        m_nphotpul[i] = v;
        totnpht += v;
    }
    
    if (totnpht != nphoton)
    {
        std::vector<long> sorter(npuls);
        
        //JCT AGAIN a patch I don't understand why "long noff = " 
        //creates a problem the line after
        int noff = nphoton - totnpht;
        int  sign = noff / abs(noff);
        
        GRBobsUtilities::getSorter(m_nphotpul, sorter);
        std::reverse(sorter.begin(), sorter.end());
        noff = abs(noff);
        
        for (long ioff=0; ioff<noff; ++ioff)
            m_nphotpul[sorter[ioff]] += sign;
    }
}




// getTmax(HepRandomEngine *engine, const int npuls, const double duration)
// Returns a vector of random numbers multiplied by the duration sorted in ascending order.  
// Random numbers are in the range [0,1).  This vector is used in the calculations of GRB times.
//
// Input:
//		engine						:	pointer to a HepRandomEngine object
//		npuls						:	number of pulses in the current burst
//		duration					:	current burst duration
//
// Output:
//		m_tmax						:	vector of (randomNumber * duration) sorted in ascending order
//
// Calls:
//		GRBobsUtilities::randGen
//		GRBobsUtilities::multiplier
//
// Caller:
//		data

void GRBpulse::getTmax(HepRandomEngine *engine, const int npuls, const double duration)
{
    m_tmax.resize(npuls);
    std::transform(m_tmax.begin(), m_tmax.end(), m_tmax.begin(), GRBobsUtilities::randGen(engine));
    
    std::transform(m_tmax.begin(), m_tmax.end(), m_tmax.begin(), GRBobsUtilities::multiplier(duration));
    
    std::sort(m_tmax.begin(), m_tmax.end());
}




// getPulse(const int npuls)
// Returns a vector needed in the calculations of GRB times.
//
// Input:
//		grbcst::nu					:	long
//		npuls						:	number of pulses in the current burst
//		m_tdiff						:	vector returned by createSigmaTdiff method
//		m_sigma						:	vector returned by createSigmaTdiff method
//
// Output:
//		m_pulse						:	output vector
//
// Calls:
//		GRBobsUtilities::multiplier
//
// Caller:
//		data

void GRBpulse::getPulse(const int npuls)
{
    DoubleSize sz = m_tdiff.size();
    
    std::vector<double> temp(sz);
    
    for (int i=0; i<sz; ++i)
        temp[i] = exp(-pow((abs(m_tdiff[i]) / m_sigma[i]), grbcst::nu));
    
    m_pulse.reserve(npuls);
    for (int ipuls=0; ipuls<npuls; ++ipuls)
    {
        std::transform(temp.begin(), temp.end(), temp.begin(), GRBobsUtilities::multiplier(m_amplitude[ipuls]));
        m_pulse.push_back(temp);
    }
}




// data(HepRandomEngine *engine, const double ethres, const long nphoton, const int npuls, const double duration)
// Returns pulse data needed in the
//
// Input:
//		engine						:	pointer to a HepRandomEngine object
//		nphoton						:	number of photons in the current burst
//		npuls						:	number of pulses in the current burst
//		duration					:	current burst duration
//
// Output:
//		deltbinsleft				:	value used in the calculations of GRB times
//		m_tmax						:	vector of (randomNumber * duration) sorted in ascending order
//		m_pulse						:	
//		m_nphotpul					:	vector containing number of photons in each pulse
//		m_sigma						:	vector containing values equal to sigrise and sigdecay
//		m_tdiff						:	
//
// Calls:
//		getAmplitude
//		getTmax
//		pickWidth
//		createSigmaTdiff
//		getPulse
//		getNphotpul
//
// Caller:
//		GRBmaker::makeTimes

long GRBpulse::data(HepRandomEngine *engine, const double ethres, const long nphoton, const int npuls, const double duration)
{
    getAmplitude(engine, npuls);
    getTmax(engine, npuls, duration);
    
    pickWidth(engine, ethres, duration);
    long deltbinsleft = createSigmaTdiff(engine);
    getPulse(npuls);
    
    getNphotpul(nphoton, npuls);
    
    return deltbinsleft;
}
