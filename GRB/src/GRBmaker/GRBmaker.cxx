// FILE: GRBmaker.cxx

// This class is instantiated in three ways:
// -  No input
//       In this mode, it generates n number of bursts and for each burst, it creates and records a photon list (time,energy)
// -  Input: Filename
//		 In this mode, it will read the photon list (time,energy) generated by the first option
// -  Input: duration, flux, fraction, power law index, npulse, flag
//		 In this mode, it creates a photon list for the burst specified by the input parameters and if the flag is set,
//			records it in a file
//

#include "GRBmaker.h"
#include "GRBpulse.h"
#include <fstream>
#include <iomanip>
#include "GRBobsConstants.h"
#include "GRBsimvecCreator.h"
#include <algorithm>              // for transform
#include <numeric>                // for accumulate
#include "CLHEP/Random/RandFlat.h"

using namespace grbobstypes;




// Default constructor  GRBmaker()
// Generate "nbsim" number of bursts.  Loops through each burst and create a photon list (time,energy) and record it in an output 
//		file.
//
// Input:
//		grbcst::seed				:	long
//		grbcst::nbsim				:	long
//
// Output:
//		GRBmaker object
//
// Calls:
//		HepRandom::getTheEngine		:	returns an engine such that engine->flat() generates random numbers between [0,1)
//		HepRandom::setTheSeed		:	sets the initial seed for the engine
//		baseFilename				:	returns the base to be used to generate the name of the output file
//		direction					:	returns the burst direction
//		grbGlobal					:	calls procedures which return samples from global distributions for GRBs: 
//											duration, peak fluxes, and power-law indices
//		makeGRB						:	computes number of photons for the current burst
//		operator<<					:	records the photon list in specified output file
//		outputFilename				:	uses base created by baseFilename to generate the name of the output file
//		specnorm					:   calculates spectral normalization
//		nphoton						:	uses specnorm and duration to calculate number of photons in the current burst

GRBmaker::GRBmaker()
	:m_univFWHM(0.0),
	 m_duration(0.0),
	 m_grbdir(std::make_pair<float,float>(0.0,0.0)),
	 m_flux(0.0),
	 m_fraction(0.0),
	 m_powerLawIndex(0.0),
	 m_specnorm(0.0),
	 m_npuls(0),
 	 m_nphoton(0),
	 m_time(),
	 m_energy()
{
	// Get an engine and set its seed
	HepRandomEngine *engine = HepRandom::getTheEngine();
	HepRandom::setTheSeed(grbcst::seed);

	// Calculate duration for each of the burst to be manufactured
	std::vector<double> duration(grbcst::nbsim);
	std::vector<double> flux(grbcst::nbsim);
	std::vector<double> fraction(grbcst::nbsim);
	std::vector<double> powerLawIndex(grbcst::nbsim);

	grbGlobal(engine, duration, flux, fraction, powerLawIndex);
	
	// create base name for the output files to be generated
	std::string base = baseFilename();
	std::string::size_type baseSize = base.size();


	// For each burst, create its direction, spectral normalization, photon energies and corresponding times
	for (long isim=0; isim<grbcst::nbsim; ++isim)
	{
		std::cout << "Processing Burst: " << isim << std::endl;

		m_specnorm = specnorm(engine, fraction[isim], powerLawIndex[isim]);

		if ((m_nphoton = nphoton(duration[isim])) > 5)
		{
			// Create burst direction
			m_grbdir = direction(engine);

			// Create spectral normalization, photon energies and corresponding times for current burst
			m_duration      = duration[isim];
			m_flux          = flux[isim];
			m_fraction      = fraction[isim];
			m_powerLawIndex = powerLawIndex[isim];
			m_npuls         = 0;
		
			makeGRB(engine);

			std::string fname = outputFilename(base, isim);
			std::ofstream os(fname.c_str());
			os << *this;
		}

		else
		{
			std::cout << "no GRB made: Insufficient flux;" << std::endl;
			std::cout << "nphoton: " << m_nphoton << std::endl << std::endl;
		}
	}
}




// Constructor  GRBmaker(filename)
// Reads the photon list (time,energy) generated by an instantiation of GRBmaker().
//
// Input:
//		filename					:	string
//
// Output:
//		GRBmaker object
//
// Calls:
//		operator>>					:	reads photon list from file specified by filename

GRBmaker::GRBmaker(const std::string &filename)
	:m_univFWHM(0.0),
	 m_duration(0.0),
	 m_grbdir(std::make_pair<float,float>(0.0,0.0)),
	 m_flux(0.0),
	 m_fraction(0.0),
	 m_powerLawIndex(0.0),
	 m_specnorm(0.0),
	 m_npuls(0),
 	 m_nphoton(0),
	 m_time(),
	 m_energy()
{
	try
	{
		std::cout << "Reading photon list from " << filename << "..." << std::endl;
		std::ifstream is (filename.c_str());
		is >> *this;

		std::ofstream os("GRB_c2.lis");
		os << *this;
	}

	catch (...)
	{
		std::cout << "Error while opening input file: " << filename << std::endl;
	}
}




// Constructor  GRBmaker(duration, npuls, flux, fraction, powerLawIndex, flag)
// Creates the photon list (time,energy) for the burst specified by the input parameters and records it in a file if the flag
//		is set.
//
// Input:
//		grbcst::seed				:	long
//
// Output:
//		GRBmaker object
//
// Calls:
//		HepRandom::getTheEngine		:	returns an engine such that engine->flat() generates random numbers between [0,1)
//		HepRandom::setTheSeed		:	sets the initial seed for the engine
//		direction					:	returns the burst direction
//		makeGRB						:	computes number of photons for the current burst
//		specnorm					:   calculates spectral normalization
//		nphoton						:	uses specnorm and duration to calculate number of photons in the current burst
//		operator<<					:	records the photon list in specified output file

GRBmaker::GRBmaker(double duration, int npuls, double flux, double fraction, double powerLawIndex, bool flag)
	:m_univFWHM(0.0),
	 m_duration(duration),
	 m_grbdir(std::make_pair<float,float>(0.0,0.0)),
	 m_flux(flux),
	 m_fraction(fraction),
	 m_powerLawIndex(powerLawIndex),
	 m_specnorm(0.0),
	 m_npuls(npuls),
 	 m_nphoton(0),
	 m_time(),
	 m_energy()
{
	HepRandomEngine *engine = HepRandom::getTheEngine();
	HepRandom::setTheSeed(grbcst::seed);

	m_specnorm = specnorm(engine, fraction, powerLawIndex);

	if ((m_nphoton = nphoton(duration)) > 5)
	{
		m_grbdir = direction(engine);

		std::cout << "Creating photon list..." << std::endl;
		makeGRB(engine);

		if (flag)
		{
			std::string out("GRB_c3.lis");
			std::ofstream os(out.c_str());
			os << *this;
		}
	}

	else 
	{
		std::cout << "no GRB made: Insufficient flux;" << std::endl;
		std::cout << "nphoton: " << m_nphoton << std::endl;
	}
}




// clone()
// Makes a copy of itself and returns it to the caller.
//
// Input:
//		this						:	pointer to the GRBmaker object to be duplicated
//
// Output:
//		GRBmaker object				:	pointer to the newly created copy of the object pointed to by the "this" pointer
//
// Calls:
//		GRBmaker(const GRBmaker &)	:	default generated copy constructor
//
// Caller:
//     interface to the simulation

GRBmaker *GRBmaker::clone() const
{
	return new GRBmaker(*this);
}




// specnorm(engine, fraction, powerLawIndex)
// Generates samples from global distribution for GRBs: duration, peak fluxes and power-law indices.
//
// Input:
//		engine						:	pointer to a HepRandomEngine object
//		fraction					:	vector of ratios of peak fluxes to maximum for "nbsim" bursts
//		powerLawIndex				:	vector containing power law indices for "nbsim" bursts
//		grbcst::ethres				:	double
//		grbcst::logdyn				:	double
//
// Output:
//		specnorm					:	spectral normalization 
//
// Calls:
//      --
//
// Caller:
//		Constructors

double GRBmaker::specnorm(HepRandomEngine *engine, double fraction, double powerLawIndex)
{
	double specnorm;

	specnorm = 3.0 * pow((grbcst::ethres*1000.), -1) * pow(fraction, 1.5);
	specnorm *= (25.0 / exp(powerLawIndex));
	specnorm *= (pow(10., (- grbcst::logdyn * engine->flat())));

	return specnorm;
}




// nphoton(duration)
// Returns number of photons in the current burst
//
// Input:
//		m_specnorm					:	spectral normalization 
//		duration					:	duration for the current burst
//
// Output:
//		m_nphoton					:	number of photons generated by the current burst
//
// Calls:
//		---
//
// Caller:
//		Constructors

long GRBmaker::nphoton(double duration)
{
	return long(m_specnorm * duration * 282743/7.0);
}




// grbGlobal(engine, duration, flux, fraction, powerLawIndex)
// Generates samples from global distribution for GRBs: duration, peak fluxes and power-law indices.
//
// Input:
//		engine						:	pointer to a HepRandomEngine object
//		grbcst::nbsim		 		:	long
//
// Output:
//		duration					:	vector containing duration for "nbsim" bursts
//		flux						:	vector containing peak fluxes for "nbsim" bursts
//		fraction					:	vector of ratios of peak fluxes to maximum for "nbsim" bursts
//		powerLawIndex				:	vector containing power law indices for "nbsim" bursts
//
// Calls:
//		getDuration					:	returns a vector of durations (size "nbsim")
//		getFlux						:	returns a vector of peak fluxes (size "nbsim")
//		getPowerLawIndex			:	returns a vector of power law indices (size "nbsim")
//
// Caller:
//		GRBmaker()					:	default constructor

void GRBmaker::grbGlobal(HepRandomEngine *engine, std::vector<double> &duration, std::vector<double> &flux, 
						 std::vector<double> &fraction, std::vector<double> &powerLawIndex)
{
	long nlong = 0L;
	getDuration(engine, nlong, duration);

	getFlux(engine, nlong, flux, fraction);

	getPowerLawIndex(engine, powerLawIndex);

	// free up memory allocated to install data for the calculations of duration, flux and power law indices
	GRBsimvecCreator::kill();	
}




// Default constructor  baseFilename()
// Returns a string containing base name used for the generation of output file name.
//
// Input:
//		grbcst::nbsim				:	long
//
// Output:
//		base						:	string used to generate the output file name
//
// Calls:
//		--
//
// Caller:
//		GRBmaker()					:	default constructor

std::string GRBmaker::baseFilename() const
{
	char *ind = new char(80);

	sprintf(ind, "%ld", grbcst::nbsim);
	int i = strlen(ind);

	std::string base = "GRB_";
	for (int j=0; j<i; ++j)
		base += '0';

	return base;
}




// Default constructor  outputFilename()
// Uses base string generated by baseFilename to create name of the output file.
//
// Input:
//		base						:	base name
//		isim						:	index for the newly generated burst
//
// Output:
//		fname						:	name of the output file
//
// Calls:
//		--
//
// Caller:
//		GRBmaker()					:	default constructor

std::string GRBmaker::outputFilename(const std::string &base, const long isim) const
{
	static std::string::size_type baseSize = base.size();

	std::string fname = base;

	char *ind = new char(80);
	sprintf(ind, "%ld", isim);

	fname.replace(baseSize-strlen(ind), baseSize-1, ind);
	fname += ".lis";

	return fname;
}




// index(engine, diff, minval, in)
// returns the index i to the last element of "in" vector such that in[i] < some random value.
//
// Input:
//		engine						:	pointer to a HepRandomEngine object
//		in							:	vector to be scanned for some specified value - SORTED IN ASCENDING ORDER
//		diff						:	difference between maximum and minimum elements of vector "in"
//		minval						:	minimum element of vector "in"
//
// Output:
//		index						:	index such that in[index] < some random value
//
// Calls:
//		engine->flat				:	returns a random number between [0,1)
//
// Caller:
//		computeFlux
//		evaluate

long GRBmaker::index(HepRandomEngine *engine, const long diff, const long minval, const std::vector<long> &in) const
{
	bool found=0;
	long i;

	while (!found)
	{
		// use random number to generate a threshold value
		const long value = long(engine->flat() * diff + minval);

		// pick last index i into "in" such that in[i] < value
		for (i=in.size()-1; i>=0 && !found; --i)
			if (in[i] < value)
				found = 1;
	}

	return i+1;
}




// evaluate(engine, diff, minval, in, v)
// Picks hi and lo values from the vector "v" and interpolates these to generate the return value.
//
// Input:
//		engine						:	pointer to a HepRandomEngine object
//		v							:	hi and lo values picked from vector "v" are interpolated to get the return value
//		in							:	vector to aid in picking values from the "v" vector
//		diff						:	difference between maximum and minimum elements of vector "in"
//		minval						:	minimum element of vector "in"
//
// Output:
//		value						:	interpolated value generated from hi and lo picked from vector "v"
//
// Calls:
//		engine->flat				:	returns a random number between [0,1)
//		GRBobsUtilities::index		:	returns a pointer into the vector "v" pointing to the lo value
//
// Caller:
//	    getDuration
//		getPowerLawIndex

double GRBmaker::evaluate(HepRandomEngine *engine, const long diff, const long minval, const std::vector<long> &in,
						  const std::vector<double> &v) const
{
	// find index loIndex such that in[loIndex] < some random number
	long loIndex = index(engine, diff, minval, in);

	double value_lo = v[loIndex+1];
	double value_hi = v[loIndex+2];

	return value_lo + engine->flat()*(value_hi - value_lo);
}




// direction(engine)
// Calculates the direction of the burst.
//
// Input:
//		engine							:	pointer to a HepRandomEngine object
//
// Output:
//		direction <zenith,azimuth> pair	:	direction of tbe current burst
//
// Calls:
//		engine->flat					:	returns a random number between [0,1)
//
// Caller:
//		constructors

std::pair<float,float> GRBmaker::direction(HepRandomEngine *engine) const
{
	float cos_zenith = 1.0 - engine->flat() * grbcst::zenNorm;

	float zenith     = acos(cos_zenith) * 180. / M_PI;
	float azimuth    = 2. * M_PI * engine->flat() * 180. / M_PI;

	return std::make_pair<float,float> (cos_zenith, azimuth);
}



	 
// computeFlux(engine, diff, minval, in, v)
// Picks hi and lo values from the vector "v" and uses the method "result" to interpolates these to generate the return value.
//
// Input:
//		engine						:	pointer to a HepRandomEngine object
//		v							:	hi and lo values picked from vector "v" are interpolated to get the return value
//		in							:	vector to aid in picking values from the "v" vector
//		diff						:	difference between maximum and minimum elements of vector "in"
//		minval						:	minimum element of vector "in"
//
// Output:
//		value						:	interpolated value generated from hi and lo picked from vector "v"
//
// Calls:
//		GRBobsUtilities::index		:	returns a pointer into the vector "v" pointing to the lo value
//		GRBobsUtilities::result		:	interpolates the lo and hi values picked from vector "v"
//
// Caller:
//		getFlux
//		

double GRBmaker::computeFlux(HepRandomEngine *engine, const long diff, const long minval, const std::vector<long> &in,
							 const std::vector<double> &v) const
{
	// find index loIndex such that in[loIndex] < some random number
	long loIndex = index(engine, diff, minval, in);

	double flux_hi = v[loIndex];
	double flux_lo = v[loIndex+1];

	return GRBobsUtilities::result(engine, flux_lo, flux_hi, grbcst::alpha);
}




// getDuration(engine, nlong, duration)
// Chooses durations from the BATSE bimodal duration distribution, where the measurement process is described by 
//		Bonnell et al. (1997, ApJ, 490, 79).  The parent sample is same as for peak fluxes: from GRB 910421 (trig #105)
//		to GRB 990123 (trig #7343).  This partial sample (1262) includes bursts where backgrounds could be fitted, and
//		peak fluxes subsequently measured.  The sample spans 7.75 years.
//
// Input:
//		grbcst::nbsim				:	long
//		engine						:	pointer to a HepRandomEngine object
//
// Output:
//		nlong						:	number of long bursts (with duration > 2.0)
//		duration					:	vector of durations for "nbsim" bursts
//
// Calls:
//		cumulativeSum				:	returns cumulative sum of the input array
//		GRBsimvecCreator::instance	:	singleton - loads data needed to generate duration, fluxes and power law indices
//		evaluate					:	evaluates the duration for current burst
//
// Caller:
//		grbGlobal
//		

void GRBmaker::getDuration(HepRandomEngine *engine, long &nlong, std::vector<double> &duration) const
{
	std::vector<double>  loEdges = GRBsimvecCreator::instance()->dur_loEdge();	
	std::vector<int>     longs   = GRBsimvecCreator::instance()->dur_long();
	std::vector<int>     shorts  = GRBsimvecCreator::instance()->dur_short();

	// Create a sum of longs and shorts in a new vector ndurs
	std::vector<int> ndurs;
	std::transform(longs.begin(), longs.end(), shorts.begin(), std::back_inserter(ndurs), std::plus<int>());

	std::vector<long> intgdurdist;
	GRBobsUtilities::cumulativeSum(ndurs, intgdurdist);

	// because intgdurdist is a cumulative sum of histplaw elements, we know that its 0th element is the minimum and
	//		last element is the max
	long diff = intgdurdist[intgdurdist.size()-1] - intgdurdist[0];
	double dur=0;

	for (int isim=0; isim<grbcst::nbsim; ++isim)
	{
		duration[isim] = (dur=evaluate(engine, diff, intgdurdist[0], intgdurdist, loEdges));
		if (dur > 2.0)
			++nlong;
	}
}




// getFLux(engine, nlong, flux, fraction)
// Chooses peak fluxes from the BATSE log N - log P;  see Bonnell et al. 1997, ApJ, 490, 79, which durplicates the procedure
//		specified by Pendleton*.  The measurement procedure is applied uniformly for that part of the BATSE sample from
//		GRB 910421 (trig #105) to GRB 990123 (trig #7343).  
//		(*Pendleton used a different PF estimation technique for the initial BATSE Catalog).
//
// Input:
//		grbcst::nbsim				:	long
//		engine						:	pointer to a HepRandomEngine object
//		nlong						:	number of long bursts (with duration > 2.0)
//											used to determine which peak flux distribution to choose from, {N,P} for longs,
//											{M,Q} for shorts
//
// Output:
//		flux						:	vector of peak fluxes for "nbsim" bursts
//		fraction					:	vector of peak flux ratios, normalized to the brightest burst in one year for "nbsim"
//											bursts
//
// Calls:
//		GRBsimvecCreator::instance	:	singleton - loads data needed to generate duration, fluxes and power law indices
//		computeFlux					:	evaluates the flux for current burst
//
// Caller:
//		grbGlobal
//		

void GRBmaker::getFlux(HepRandomEngine *engine, long nlong, std::vector<double> &flux, 
					   std::vector<double> &fraction) const
{
	std::vector<long>    m = GRBsimvecCreator::instance()->flux_m();
	std::vector<long>    n = GRBsimvecCreator::instance()->flux_n();
	std::vector<double>  p = GRBsimvecCreator::instance()->flux_p();
	std::vector<double>  q = GRBsimvecCreator::instance()->flux_q();

	grbobstypes::DoubleConstIter it0 = std::max_element(p.begin(), p.end());
	grbobstypes::DoubleConstIter it1 = std::max_element(q.begin(), q.end());
	double maxP = std::max<double> (*it0, *it1);

	grbobstypes::LongConstIter it2 = std::max_element(n.begin(), n.end());
	grbobstypes::LongConstIter it3 = std::min_element(n.begin(), n.end());
	long diff = *it2 - *it3;

	double f=0;

	for (long isim=0; isim<nlong; ++isim)
	{
		flux[isim]     = (f=computeFlux(engine, diff, *it3, n, p));
		fraction[isim] = (f/maxP);
	}

	it2 = std::max_element(m.begin(), m.end());
	it3 = std::min_element(m.begin(), m.end());
	diff = *it2 - *it3;

	for (isim=nlong; isim<grbcst::nbsim; ++isim)
	{
		flux[isim]     = (f=computeFlux(engine, diff, *it3, m, q));
		fraction[isim] = (f/maxP);
	}
}




// getPowerLawIndex(engine, powerLawIndex)
// Chooses spectral power-law indices from the BATSE power-law distribution, as measured by Preece et al. (1999)
//
// Input:
//		grbcst::nbsim				:	long
//		engine						:	pointer to a HepRandomEngine object
//
// Output:
//		powerLawIndex				:	vector of power law indices for "nbsim" bursts
//
// Calls:
//		cumulativeSum				:	returns cumulative sum of the input array
//		GRBsimvecCreator::instance	:	singleton - loads data needed to generate duration, fluxes and power law indices
//		evaluate					:	evaluates the power law index for current burst
//
// Caller:
//		grbGlobal
//		

void GRBmaker::getPowerLawIndex(HepRandomEngine *engine, std::vector<double> &powerLawIndex) const
{
	std::vector<double>  loEdges  = GRBsimvecCreator::instance()->pl_loEdge();
	std::vector<int>     histplaw = GRBsimvecCreator::instance()->pl_histplaw();

	std::vector<long>  intgplawdist;
	GRBobsUtilities::cumulativeSum(histplaw, intgplawdist);

	// because intgplawdist is a cumulative sum of histplaw elements, we know that its 0th element is the minimum and
	//		last element is the max
	long diff = intgplawdist[intgplawdist.size()-1] - intgplawdist[0];

	for (long isim=0; isim<grbcst::nbsim; ++isim)
		powerLawIndex[isim] = (evaluate(engine, diff, intgplawdist[0], intgplawdist, loEdges));
}




// getTimes(engine, nphots, deltbinsleft, iphotoff, tmax, pulse)
// Uses pulse data to generate list of photon times
//
// Input:
//		grbcst::timres				:	double
//		grbcst::ethres				:	double
//		engine						:	pointer to a HepRandomEngine object
//		nphots						:	number of photons in the current pulse
//		deltbinsleft				:	value generated by GRBpulse::createSigmaTdiff method
//		iphotoff					:	index to the first photon to generate time for
//		tmax						:	time for peak amplitude for the current burst and pulse
//		pulse						:	vector generated by GRBpulse::getPulse method
//
// Output:
//		m_times						:	"nphots" photon times for the elements in range [iphotoff,iphotoff+nphots)
//
// Calls:
//		cumulativeSum				:	returns cumulative sum of the input array
//		engine->flat				:	returns a random number between [0,1)
//
// Caller:
//		makeTimes
//		

void GRBmaker::getTimes(HepRandomEngine *engine, const long nphots, const long deltbinsleft, const long iphotoff, 
						const double tmax, const std::vector<double> &pulse)
{
	std::vector<double> cumpulse;
	
	GRBobsUtilities::cumulativeSum(pulse, cumpulse);
	double maxval = cumpulse[cumpulse.size()-1];
	std::transform(cumpulse.begin(), cumpulse.end(), cumpulse.begin(), GRBobsUtilities::multiplier(1./maxval));

	double photontime, efactor, twidthscale;

	for (int iphot=0; iphot<nphots; ++iphot)
	{
		grbobstypes::DoubleIter it = cumpulse.end();
		while (it == cumpulse.end())
			it = std::find_if(cumpulse.begin(), cumpulse.end(), std::bind2nd(std::greater_equal<double>(), engine->flat()));

		int index = std::distance(cumpulse.begin(), it);

		photontime   = (index - deltbinsleft) * grbcst::timres;
		efactor      = pow((m_energy[iphot+iphotoff] / grbcst::ethres), -0.333);
		twidthscale  = photontime * efactor;

		m_time[iphot+iphotoff] = twidthscale + tmax + 0.5 * efactor * deltbinsleft * grbcst::timres;
	}
}




// makeTimes(engine)
// Makes BATSE-like GRB time profiles, placing GLAST photons a la cumulative BATSE intensity, but in narrower pulses:
//		(1) m_npuls = number of pulses, proportional to BATSE duration.
//		(2) pulse peak amplitude is random (0.0=>1.0); sort amps in descending amp order.
//		(3)	scramble amps of {1st,2nd} halves of pulses, separately (leaves profile assymmetric).
//		(4)	center of pulse time is random within duration.  sort the times in ascending order.
//		(5)	pulse width is drawn from BATSE width distribution for bright bursts (attributes paper), scaled to GLAST energies.
//				using width ~E^-0.333.
//		(6)	make m_npuls pulses with "bisigma" shapes => sum to produce time profile
//		(7)	form cumulative distribution of BATSE intensity
//		(8) distribute the m_nphoton photons according to cumulative intensity => GRBtimes
//		(9)	offset the photon times according to (a) energy dependence, width ~E^-0.333 and (b) time of peak, also 
//				proportional to E^-0.333.
//
// Input:
//		engine						:	pointer to a HepRandomEngine object
//		m_npuls						:	0:  generate number of pulses
//										>0: use this number to generate pulse data
//
// Output:
//		m_npuls						:	number of pulses
//		m_times						:	photon times sorted in ascending order
//
// Calls:
//		GRBpulse::data				:	returns pulse data for current burst -
//											- tmax		:	times of peak amplitudes
//											- pulse		:	vector created by GRBpulse::getPulse
//											- nphotpul	:	number of photons in each pulse
//		getTimes					:	returns photon times for nphotpul[ipulse] where 0<=ipulse<m_npuls
//		GRBpulse::univFWHM			:	returns the universal width
//
// Caller:
//		makeGRB

void GRBmaker::makeTimes(HepRandomEngine *engine)
{
	if (m_npuls == 0)
		m_npuls = std::max(int(15 * m_duration/30 + 0.5), 1);

	try
	{
		GRBpulse *grbPulse = new GRBpulse;

		long deltbinsleft = grbPulse->data(engine, m_nphoton, m_npuls, m_duration);

		// create times for each pulse
		const std::vector<double> tmax                 = grbPulse->tmax();
		const std::vector< std::vector<double> > pulse = grbPulse->pulse();
		const std::vector<long> nphotpul               = grbPulse->nphotpul();

		long iphotoff = 0;
		for (long ipuls=0; ipuls<m_npuls; ++ipuls)
		{
			if (nphotpul[ipuls] > 0)
			{
				getTimes(engine, nphotpul[ipuls], deltbinsleft, iphotoff, tmax[ipuls], pulse[ipuls]);
				iphotoff += nphotpul[ipuls];
			}
		}
		
		m_univFWHM = grbPulse->univFWHM();

		delete grbPulse;
	}

	catch(...)
	{
		std::cout << "Cannot allocate memory for the GRB pulse..." << std::endl;
	}
}




// makeEnergies(engine)
// Computes photon energies for the current burst: 
//
// Input:
//		engine						:	pointer to a HepRandomEngine object
//		grbcst::ethres				:	double
//		grbcst::logdyn				:	double
//		grbcst::emax				:	double
//		m_powerLawIndex				:	power law index for the current burst
//
// Output:
//		m_energy					:	photon energies sorted in time order
//
// Calls:
//		engine->flat				:	returns a random number between [0,1)
//
// Caller:
//		makeGRB

void GRBmaker::makeEnergies(HepRandomEngine *engine)
{
	if (m_powerLawIndex == 1.0)
	{
		for (long iphot=0; iphot!=m_nphoton; ++iphot)
			m_energy[iphot] = grbcst::ethres * exp(engine->flat() * log(grbcst::emax/grbcst::ethres)); 
	}

	else
	{
		double x = 1.0 - exp((1.0 - m_powerLawIndex) * log(grbcst::emax/grbcst::ethres));

		for (long iphot=0; iphot!=m_nphoton; ++iphot)
			m_energy[iphot] = grbcst::ethres * exp(log(1.0 - x*engine->flat()) / (1.0 - m_powerLawIndex));
	}
}




// makeGRB(engine)
// Computes number of photons for the current burst:
//		SpecNorm units, integrated Peak Flux: photons cm^-2 s^-1 (>grbcst::ethres).
//		Researches on normalization:
//			(1) Bonnell's fits to bright BATSE bursts.
//			(2)	comparison with EGRET norms for bright bursts - Catelli's, Dingus' and Schneid's works.
//			(3) analysis of Preece et al. spectroscopy catalog of bright BATSE bursts (see JPN routine Specanal.pro).
//
//		The cofactors for SpecNorm are:
//			(a) {average flux/peak flux} = ~1/7.
//			(b)	scaling by (peak_flux)^1.5, determined from inspection of Preece et al..
//			(c)	duration (seconds).
//			(d)	282743 cm^2 (6-meter dia. illuminated disk).
//			(e)	scaling to integral above Ethres (e.g., 0.03 GeV) for case beta=-2.
//			(f)	dispersion (dynrange) to approximately replicate the scatter in peak flux vs. normalization at 1 MeV
//					as estimated from Preece et al. catalog.
//			(g)	a dependence on power-law index as estimated from Preece et al. catalog.
//
//		Thus, m_nphoton is number of photons normally incident on projected disk of GLAST illumination sphere, integrated
//			above Ethres, for chosen peak flux and duration.  Energies distributed as power-law index.
// 
//
// Input:
//		engine						:	pointer to a HepRandomEngine object
//
// Output:
//		m_times						:	photon times sorted in ascending order
//		m_energies					:	photon energies sorted in time order
//
// Calls:
//		engine->flat				:	returns a random number between [0,1)
//		makeEnergy					:	makes photon energies
//		makeTimes					:	makes BATSE-like GRB time profiles
//		GRBobsUtilities::sortVector	:	sorts vector A in the order of vector B
//
// Caller:
//		makeTimes

void GRBmaker::makeGRB(HepRandomEngine *engine)
{
	m_energy.resize(m_nphoton);
	m_time.resize(m_nphoton);

	makeEnergies(engine);

	makeTimes(engine);

	// sort GRB times
	std::vector<double> time(m_time);
	std::sort(m_time.begin(), m_time.end());

	// sort GRB energies in time order
	GRBobsUtilities::sortVector(0, time, m_time, m_energy);
}




// operator>>(is, grbMaker)
// Reads burst information and photon list (time,energy) from the input stream and updates the object grbMaker with it
//
// Input:
//		is							:	input stream
//
// Output:
//		is							:	input stream
//		grbMaker					:	GRBmaker object containing following information:
//			m_nphoton					:	number of photons generated by the current burst
//			m_flux						:	peak flux
//			m_duration					:	burst duration 
//			m_powerLawIndex				:	power law index
//			m_specnorm					:	spectral normalization 
//			m_npuls						:	number of pulses in the current burst
//			m_times						:	photon times 
//			m_energy					:	photon energies 
//
// Calls:
//		---
//
// Caller:
//		constructor GRBmaker(filename)

std::ifstream &operator>>(std::ifstream &is, GRBmaker &grbMaker)
{
	std::string str;
	is >> str >> str >> str >> str >> str >> str;
	is >> grbMaker.m_flux >> grbMaker.m_duration >> grbMaker.m_powerLawIndex >> grbMaker.m_specnorm >> grbMaker.m_npuls >>
		grbMaker.m_univFWHM;
	is >> str >> str >> str;

	float zenith, azimuth;
	is >> zenith >> azimuth;
	grbMaker.m_grbdir = std::make_pair<float,float> (zenith, azimuth);

	is >> str >> str;
	is >> grbMaker.m_nphoton;
	is >> str >> str >> str >> str;

	grbMaker.m_time.resize(grbMaker.m_nphoton);
	grbMaker.m_energy.resize(grbMaker.m_nphoton);

	for (long i=0; i<grbMaker.m_nphoton; ++i)
		is >> grbMaker.m_time[i] >> grbMaker.m_energy[i];

	return is;
}




// operator>>(os, grbMaker)
// Writes burst information and photon list (time,energy) contained in the object grbMaker to the output stream
//
// Input:
//		is							:	input stream
//		grbMaker					:	GRBmaker object containing following information:
//			m_nphoton					:	number of photons generated by the current burst
//			m_flux						:	peak flux
//			m_duration					:	burst duration 
//			m_powerLawIndex				:	power law index
//			m_specnorm					:	spectral normalization 
//			m_npuls						:	number of pulses in the current burst
//			m_times						:	photon times 
//			m_energy					:	photon energies 
//
// Output:
//		is							:	input stream
//
// Calls:
//		---
//
// Caller:
//		default constructor GRBmaker()

std::ofstream &operator<<(std::ofstream &os, const GRBmaker &grbMaker)
{
	std::pair<double,double> grbdir = grbMaker.dir();

	os << "Fp, dur, beta, Specnorm, Npulses, UnivFWHM= " << std::setw(12) << std::fixed << grbMaker.flux() << "   " <<
		grbMaker.duration() << "   " << grbMaker.powerLawIndex() << "   " << grbMaker.specnorm() << "   " << 
		grbMaker.npuls() << "  " << grbMaker.univFWHM() << std::endl;

	os << "ZenAng, AziAng = " << std::setw(12) << std::fixed << grbdir.first << "   " << grbdir.second << std::endl;

	os << "nphotons = " << grbMaker.nphoton() << std::endl;

	os << "(Times, Energies) per photon:" << std::endl;


	if (grbMaker.nphoton() > 0)
	{
		std::vector<double> time   = grbMaker.time();
		std::vector<double> energy = grbMaker.energy();

		for (long i=0; i<grbMaker.nphoton(); ++i)
			os << std::setw(12) << std::fixed << time[i] << "   " << energy[i]<< std::endl;
	}

	return os;
}


