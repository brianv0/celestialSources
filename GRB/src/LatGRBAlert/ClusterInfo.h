/*!
* \class ClusterInfo
*
* \brief This class provides interface for cluster of photons.
* \author Jay Norris        jnorris@lheapop.gsfc.nasa.gov
* \author Sandhia Bansal    sandhiab@lheapop.gsfc.nasa.gov
*
*/
#ifndef CLUSTER_INFO_H
#define CLUSTER_INFO_H


class ClusterInfo
{
public:
    inline void setTime(double time)       { m_time=time; }
    inline void setEnergy(double energy)   { m_energy=energy; }
    inline void setDist(double dist)       { m_dist=dist; }
    inline void setSignal(bool signal)     { m_signal=signal; }
    
    inline double time() const     { return m_time; }
    inline double energy() const   { return m_energy; }
    inline double dist() const     { return m_dist; }
    inline bool signal() const     { return m_signal; }
    
private:
    double m_time;
    double m_energy;
    double m_dist;
    bool   m_signal;
};

#endif
