// FILE: LatGRBAlert.cxx


#include <algorithm>   // for count_if
#include <iostream>

#include "LatGRBAlert.h"
#include "PhotonInfo.h"
#include "LatGRBAlertConstants.h"
#include "ClusterData.h"

using namespace lattrigtypes;




LatGRBAlert::LatGRBAlert(const std::vector<PhotonInfo> &photonData, const long nbckoff, const long nGrb, const long nBck)
	:m_falseTrig(0),
	 m_nbckoff(nbckoff),
	 m_firstTriggerTime(0),
	 m_highbck(0),
	 m_maxResult(0),
	 m_photonData(photonData),
	 m_jointLike(),
	 m_grbTruth(),
	 m_grb1phot()
{
	long nregion = (nGrb + nBck - lattrigcst::nrange) / lattrigcst::nmove - 1;

	if (nregion <= 0)
		std::cout << "Not enough photons" << std::endl;
	
	else
	{
		m_jointLike.reserve(nregion);
		m_grbTruth.resize(nregion);
		m_grb1phot.resize(nregion);
	
		std::vector<PhotonInfo> pdata(lattrigcst::nrange);
		std::vector<long>       dumhighbck;

		ClusterData *clusterData = new ClusterData;

		long                    offset=0;

		std::cout << "number of regions: " << nregion << std::endl;
		for (long region=0; region<nregion; ++region)
		{
			std::cout << "process region " << region << std::endl;

			std::copy(&m_photonData[offset], &m_photonData[offset+lattrigcst::nrange], pdata.begin());

			long nsignal = std::count_if(pdata.begin(), pdata.end(), isSignalSet());

			m_grbTruth[region]  = ((nsignal >= 5) ? 1 : 0);
			m_grb1phot[region]  = ((nsignal >= 1) ? 1 : 0);

			m_jointLike.push_back(clusterData->triggerLikelihood(pdata));

			// Don't look for false triggers for now -
			//if (m_grb1phot[region] == 0)
			//{
			//	dumhighbck.push_back(region);
			//	if (!m_falseTrig)   
			//		m_falseTrig = isFalseTrigger(region);
			//}

			if (m_firstTriggerTime == 0)
			{
				m_firstTriggerTime = firstTriggerTime(region, pdata);

				if (m_firstTriggerTime > 0)   // found first trigger
					break;
			}

			offset += lattrigcst::nmove;
		}

		// Don't record false triggers for now -
		//recordFalseTriggers(dumhighbck);

		if (m_jointLike.empty())
			std::cout << ">>>>>> No Trigger Found. <<<<<<" << std::endl;

		delete clusterData;
	}
}



bool LatGRBAlert::isSignalSet::operator () (const PhotonInfo &x)
{
	return x.signal() == 1;
}



bool LatGRBAlert::isFalseTrigger(const long region)
{
	if ((region >= m_nbckoff) && (m_jointLike[region] > m_jointLike[region-m_nbckoff]))
		return 1;

	return 0;
}


double LatGRBAlert::firstTriggerTime(const long region, const std::vector<PhotonInfo> &photonData)
{
	if ((region >= m_nbckoff) && (m_jointLike[region] > m_jointLike[region-m_nbckoff]))
		return photonData[region].time();

	return 0;
}


void LatGRBAlert::recordFalseTriggers(const std::vector<long> &vindex)
{
	DoubleConstIter it = std::max_element(&m_jointLike[vindex[0]], &m_jointLike[vindex[vindex.size()-1]]);
	m_highbck = *it;

	it = std::max_element(m_jointLike.begin(), m_jointLike.end());
	m_maxResult = *it;
}
