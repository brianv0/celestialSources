// FILE: ClusterData.cxx


#include <algorithm>   // for sort
#include <numeric>     // for accumulate
#include <cmath>       // for log10

#include "ClusterData.h"
#include "LatTriggerConstants.h"

using namespace lattrigtypes;



#include <fstream>
static std::ofstream os2("likeinfo.lis");


// Default constructor  ClusterData()

ClusterData::ClusterData()
	: m_photonData(),
	  m_useDist(),
	  m_avgDist(),
	  m_extractedClusterData(),
	  m_jointLike(0)
{
}



// cosDistance()
// 
//
// Input:
//      index1, index2              :   indices pointing to two photons in the list
// 
// Output:
//		distance between two photons
//
// Calls:
//      ---

const double ClusterData::cosDistance(const long index1, const long index2) const
{
	double phidiff = abs(m_photonData[index1].phi() - m_photonData[index2].phi());

	return (cos(m_photonData[index1].theta()) *
			cos(m_photonData[index2].theta())) +
		   (sin(m_photonData[index1].theta()) *
		    sin(m_photonData[index2].theta()) *
		    cos(phidiff));
}



// fillUseDistance(const long index)
// form N*(N-1) distances (on sphere between nrange photons).  Use spherical trig:
// pole=zenith.  For each photon, there are (N-1) distances.  Per "photonl" sort these
// distances. we will keep the smallest "cluster" -- that distance set with the smallest
// average distance for photons within 95% containment for 30MeV photon (using PSF for layers 12-15).  
// This distance = ~2.3*sigma68=2.3*13.6 ~= 34 degrees.
//
// Input:
//      index                       :   pointer to the first photon to be considered
// 
// Output:
//		distances between photons
//
// Calls:
//      ---

void ClusterData::fillUseDistance(const long index)
{
	std::vector<double> tmpdist;
	tmpdist.reserve(lattrigcst::nrange-1);

	for (long i=index; i<lattrigcst::nrange; ++i)
	{
		tmpdist.clear();

		for (long j=0; j<lattrigcst::nrange; ++j)
		{
			if (i != j)
				tmpdist.push_back(acos(cosDistance(i, j)));
		}

		std::sort(tmpdist.begin(), tmpdist.end());
		m_useDist[i] = tmpdist;
	}
}



// calcUseDistance()
// calls fillUseDistance to form N*(N-1) distances (on sphere between nrange photons).  
// we are moving forward in "nmove" chunks, so it is only the first time that all distances need to be calculated;
// after that only following elements need to be recalculated:
//      - "nmove" elements of the first (nragne-nmove) rows
//      - last nmove rows
//
// Input:
//      lattrigcst::nrange            :   number of photons in this cluster
//      lattrigcst::nmove             :   number of photons in a chunk to move along nrange photons
//      m_useDist                     :   distances between nrange photons; empty for the first time
// 
// Output:
//      m_useDist                     :   distances between nrange photons
//
// Calls:
//      ---

void ClusterData::calcUseDistance()
{
	if (m_useDist.empty())
	{
		m_useDist.resize(lattrigcst::nrange);
		fillUseDistance(0);
	}

	else
	{
		std::copy(&m_useDist[lattrigcst::nmove], m_useDist.end(), m_useDist.begin());

		// replace first "nmove" elements of the first (nrange-nmove) rows
		long nrows = lattrigcst::nrange - lattrigcst::nmove;

		long secondIndex = m_useDist[0].size() - lattrigcst::nmove + 1;

		for (long i=0; i<nrows; ++i)
		{
			for (long j=0; j<lattrigcst::nmove; ++j)
				m_useDist[i][j] = acos(cosDistance(i, secondIndex+j));

			std::sort(m_useDist[i].begin(), m_useDist[i].end());
		}

		fillUseDistance(nrows);
	}
}



// calcAvgDistance()
// find most compact cluster, considering photons within 2 sigma of 30 MeV photon radius.
//
// Input:
//      lattrigcst::nrange            :   number of photons in this cluster
//      lattrigcst::twosigdist        :   
//      m_useDist                     :   distances between nrange photons
// 
// Output:
//      m_avgDist                     :   average distances between nrange photons
//
// Calls:
//      ---

void ClusterData::calcAvgDistance()
{
	m_avgDist.clear();
	m_avgDist.reserve(lattrigcst::nrange);

	long i=0;
	for (DoubleVecConstIter it=m_useDist.begin(); it!=m_useDist.end(); ++it)
	{
		// m_useDist[i,*] is a sorted vector
		DoubleConstRevIter iter = std::find_if((*it).rbegin(), (*it).rend(), 
			std::bind2nd(std::less<double>(), lattrigcst::twosigdist));

		long dist = std::distance(iter, (*it).rend());

		if (dist > 0)
		{
			double value = std::accumulate(iter, (*it).rend(), 0.0);
			m_avgDist.push_back(std::accumulate(iter, (*it).rend(), 0.0)/dist);
		}

		else
			m_avgDist.push_back(9999.);
	}
}



// greater_equal()
// return true if distance associated with current record is >= some specified value.
//
// Input:
//      x                             :   current cluster record
//      m_value                       :   user supplied value
// 
// Output:
//      true if x.dist >= m_value
//
// Calls:
//      ---

const bool ClusterData::greater_equal::operator () (ClusterInfo &x) const
{
	return x.dist() >= m_value;
}



// less()
// return true if distance associated with current record is < some specified value.
//
// Input:
//      x                             :   current cluster record
//      m_value                       :   user supplied value
// 
// Output:
//      true if x.dist < m_value
//
// Calls:
//      ---

const bool ClusterData::less::operator () (ClusterInfo &x) const
{
	return x.dist() < m_value;
}



// pdistValue(const double distance)
//
// Input:
//      distance                      :   input value
// 
// Output:
//      distance/(2*M_PI)
//
// Calls:
//      ---

const double ClusterData::pdistValue(const double distance) const
{
	return distance / (2 * M_PI);
}



// pdtValue(const double deltaTime)
// see JPN thesis, page 101: Prob of nu counts in interval <=dt:
// 1-exp(-x)*(1+x), where x=R*dt, and R=expected background rate=4.07 Hz (avg)
//
// Input:
//      lattrigcst::ratebck           :   average background rate
//      deltaTime                     :   input value
// 
// Output:
//      double
//
// Calls:
//      ---

const double ClusterData::pdtValue(const double deltaTime) const
{
	bool factor = ((deltaTime > 0.0005) ? 1 : 0);
	double x    = lattrigcst::ratebck * factor;

	return (1 - exp(-x) * (1+x));
}



// psfValue(const double energy, const double denom)
// energy weighting is proportional to the PSF (photon).  For now, since we didn't log fst_x_lyr values, use JTB's
// parameterization for layers 12-15: alp1=0.93 & C1=0.52 & C2=0.063.
//
// Input:
//      energy                        :   input value
//      denom                         :   input value
// 
// Output:
//      double
//
// Calls:
//      ---

const double ClusterData::psfValue(const double energy, const double denom) const
{
	return (0.52 * pow(energy, -0.93) + 0.063) / denom;
}



// likelihood(const double dist, const double dt, const double psf)
//
// Input:
//      dist, dt, psf                      :   double
// 
// Output:
//      likelihood                         :   trigger likelihood
//
// Calls:
//      ---

const double ClusterData::likelihood(const double dist, const double dt, const double psf) const
{
	os2 << "likelihood::dist: " << dist << " dt: " << dt << " psf: " << psf << "return value: " << -(dist+dt+psf) << std::endl;
	return -(dist + dt + psf);
}



// extractClusterData(const long chosenphot)
// re-extract cluster (N-1) distances, and associated photon parameters.  Fast approach appears to be: recompute
// distances from "chosen" photon with tightest cluster.
//
// Input:
//      chosenphot                         :   index to the photon with tightest cluster
// 
// Output:
//      m_extractedClusterData             :   extracted cluster data
//
// Calls:
//      greater_equal

void ClusterData::extractClusterData(const long chosenphot)
{
	m_extractedClusterData.clear();
	m_extractedClusterData.reserve(lattrigcst::nrange);
	ClusterInfo record;

	for (long i=0; i<lattrigcst::nrange; ++i)
	{
		if (i != chosenphot)
		{
			record.setDist(acos(cosDistance(chosenphot, i)));
			record.setTime(m_photonData[i].time());
			record.setEnergy(m_photonData[i].energy());
			record.setSignal(m_photonData[i].signal());

			m_extractedClusterData.push_back(record);
		}
	}

	std::sort(m_extractedClusterData.begin(), m_extractedClusterData.end(), distCmp());

	std::vector<ClusterInfo>::iterator it = 
		std::find_if(m_extractedClusterData.begin(), m_extractedClusterData.end(), 
		greater_equal(lattrigcst::twosigdist));

	m_extractedClusterData.erase(it, m_extractedClusterData.end());

	// set time, energy and signal for next reoord (not dist)
	record.setTime(m_photonData[chosenphot].time());
	record.setEnergy(m_photonData[chosenphot].energy());
	record.setSignal(m_photonData[chosenphot].signal());

	m_extractedClusterData.push_back(record);
}



// jointLikelihood()
// calculate trigger likelihood.
//
// Input:
//      ---
// 
// Output:
//      likelihood                         :   trigger likelihood
//
// Calls:
//      likelihood
//      psfValue

const double ClusterData::jointLikelihood() const
{
	std::vector<ClusterData>::size_type sz = m_extractedClusterData.size();
	std::vector<double> clusterTimes;
	clusterTimes.reserve(sz);

	for (std::vector<ClusterInfo>::const_iterator it=m_extractedClusterData.begin(); it!=m_extractedClusterData.end(); ++it)
		clusterTimes.push_back((*it).time());

	std::sort(clusterTimes.begin(), clusterTimes.end());

	// compute relative probability of chance occurrence.
	// joint probability is product of : P(dts), P(dists), P(Energies)

	double pdist = 0;
	double pdt   = 0;
	double psf   = 0;

	double denom = 0.52 * pow(0.30, -0.93) + 0.063;

	--sz;
	for (long i=0; i<sz; ++i)
	{
		pdist += (log10(pdistValue(m_extractedClusterData[i].dist())));
		pdt   += (log10(pdtValue(clusterTimes[i+1] - clusterTimes[i])));
		psf   += (log10(psfValue(m_extractedClusterData[i].energy(), denom)));
	}

	psf   += (log10(psfValue(m_extractedClusterData[i].energy(), denom)));

	return likelihood(pdist, pdt, psf);
}



// triggerLikelihood()
// calculate trigger likelihood.
//
// Input:
//      ---
// 
// Output:
//      likelihood                         :   trigger likelihood
//
// Calls:
//      calcUseDistance
//      calcUseDistance
//      extractClusterData
//      jointLikelihood

const double ClusterData::triggerLikelihood(const std::vector<PhotonInfo> &photonData)
{
	m_photonData = photonData;

	calcUseDistance();

	calcAvgDistance();

	DoubleIter it = std::min_element(m_avgDist.begin(), m_avgDist.end());

	if (*it >= 9999.)
		return 0;
	else
	{
		extractClusterData(std::distance(m_avgDist.begin(), it));
		return jointLikelihood();
	}
}