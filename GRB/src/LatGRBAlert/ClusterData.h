/*!
* \class ClusterData
*
* \brief This class provides interface for cluster of photons.
* \author Jay Norris        jnorris@lheapop.gsfc.nasa.gov
* \author Sandhia Bansal    sandhiab@lheapop.gsfc.nasa.gov
*
* This class implements methods to create distances between photons in a cluster, and
* to compute relative probability of trigger.
*
*/

#ifndef CLUSTER_DATA_H
#define CLUSTER_DATA_H

#include <vector>
#include <functional>  // for unary_function

#include "PhotonInfo.h"
#include "ClusterInfo.h"




class ClusterData
{
public:
    /*!
     * Constructor.
     */
    ClusterData();
    

    /*!
     * \brief Compute probability of a trigger.
     */
    const double triggerLikelihood(const std::vector<PhotonInfo> &photonData);
    
private:
    /*!
     * \brief Find most compact cluster, considering photons within 2 sigma of 30 MeV photon radius.
     */
    void   calcAvgDistance();

    /*! 
     * \brief Calls fillUseDistance to form N*(N-1) distances (on sphere between nrange photons). 
     */
    void   calcUseDistance();

    /*!
     * \brief Return cosine of distance.
     */
    const double cosDistance(const long index1, const long index2) const;

    /*!
     * \brief Form N*(N-1) distances (on sphere between nrange photons).  
     *
     * Use spherical trig:pole=zenith.  For each photon, there are (N-1) distances.  
     * Per "photonl" sort these distances. we will keep the smallest "cluster" -- 
     * that distance set with the smallest average distance for photons within 95% containment 
     * for 30MeV photon (using PSF for layers 12-15).  
     * This distance = ~2.3*sigma68=2.3*13.6 ~= 34 degrees.
     */
    void   fillUseDistance(const long index);

    /*!
     * \brief Extract cluster (N-1) distances, and associated photon parameters.  
     * Fast approach appears to be: recompute distances from "chosen" photon with tightest cluster.
     */
    void   extractClusterData(const long chosenphot);

    /*! 
     * \brief Compute relative probability of chance occurrence.
     */
    const double jointLikelihood() const;

    /*!
     * \breif Compute trigger likelihood from input parameters.
     */
    const double likelihood(const double dist, const double dt, 
        const double psf) const;

    /*!
     * \brief Compute one of the parameters needed to determine likelihood.
     */
    const double pdistValue(const double distance) const;

    /*!
     * \brief Compute one of the parameters needed to determine likelihood.
     */
    const double pdtValue(const double deltaTime) const;

    /*!
     * \brief Compute one of the parameters needed to determine likelihood.
     */
    const double psfValue(const double energy, const double denom) const;

    
    
    // Data members
    std::vector<PhotonInfo>           m_photonData;             //! Photon list (time,energy,theta,phi,signal)
    std::vector<std::vector<double> > m_useDist;                //! Distances between cluster points
    std::vector<double>               m_avgDist;                //! Distances needed to find most compact cluster
    std::vector<ClusterInfo>          m_extractedClusterData;   //! Extracted cluster data
    double                            m_jointLike;              //! Joint probability
    
    
    
    /*!
     * This class provides functionality to sort a vector in decreasing order of distances.
     */
    struct greater_equal : public std::unary_function<ClusterInfo, ClusterInfo>
    {
        greater_equal(double value) : m_value(value) {}
        const bool operator() (ClusterInfo &x) const;
        
        double m_value;
    };
    
    
    /*! 
     * This class provides functionality to sort a vector in increasing order of distances.
     */
    struct less : public std::unary_function<ClusterInfo, ClusterInfo>
    {
        less(double value) : m_value(value) {}
        const bool operator() (ClusterInfo &x) const;
        
        double m_value;
    };
};



/*!
 * This class implements method to sort data by distances
 */
class distCmp
{
public:
    bool operator()(ClusterInfo &data1, ClusterInfo &data2)
    {
        return data1.dist() < data2.dist();    
    }
};


#endif