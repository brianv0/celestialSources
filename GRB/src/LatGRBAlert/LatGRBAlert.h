/*!
* \class LatGRBAlert
*
* \brief This class provides interface for trigger likelihood determination process.
* \author Jay Norris        jnorris@lheapop.gsfc.nasa.gov
* \author Sandhia Bansal    sandhiab@lheapop.gsfc.nasa.gov
*
*/


#ifndef LATGRBALERT_H
#define LATGRBALERT_H

#include <vector>
#include <deque>
#include <functional>   // for unary_function

class PhotonInfo;



class LatGRBAlert
{
public:
    /*! Constructor.
     */
    LatGRBAlert(const std::vector<PhotonInfo> &photonData, const long nbckoff, 
        const long nGrb, const long nBck);
    
    

    // Accessors
    bool falseTrig() const                             { return m_falseTrig; }
    long nbckoff() const                               { return m_nbckoff; }
    double firstTriggerTime() const                    { return m_firstTriggerTime; }
    double highbck() const                             { return m_highbck; }
    double maxResult() const                           { return m_maxResult; }
    const std::vector<PhotonInfo> &photonData() const  { return m_photonData; }
    const std::vector<double> &jointLike() const       { return m_jointLike; }
    const std::deque<bool> &grbTruth() const           { return m_grbTruth; }
    const std::deque<bool> &grb1phot() const           { return m_grb1phot; }
    
    
private:
    /*!
     * \brief Checks regions past (m_nbckoff)-th region to return time of first trigger.
     */
    double firstTriggerTime(const long region, const std::vector<PhotonInfo> &photonData);

    //bool isFalseTrigger(const long region);
    //bool isSignal();
    //void recordFalseTriggers(const std::vector<long> &vindex);
    
    
    // Data members
    bool                     m_falseTrig;              //! true - false trigger - not yet used
    long                     m_nbckoff;                //! Examine regions past this number for trigger
    double                   m_firstTriggerTime;       //! Time of the first trigger
    double                   m_highbck;                //! Not yet used
    double                   m_maxResult;              //! Not yet used
    std::vector<PhotonInfo>  m_photonData;             //! Photon list
    std::vector<double>      m_jointLike;              //! List of joint likelihood
    std::deque<bool>         m_grbTruth;               //! Set to true for the region that contains >= 5 signals
    std::deque<bool>         m_grb1phot;               //! Set to true for the region that contains >= 1 signals
    
    
    
    /*!
     * \brief This class provides method to check whether the event is a signal or background
     */
    struct isSignalSet : public std::unary_function<PhotonInfo, PhotonInfo>
    {
        isSignalSet() {}
        bool operator() (const PhotonInfo &x);
    };
};


#endif