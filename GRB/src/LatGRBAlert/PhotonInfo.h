/*!
* \class PhotonInfo
*
* \brief This class provides interface for photons.
* \author Jay Norris        jnorris@lheapop.gsfc.nasa.gov
* \author Sandhia Bansal    sandhiab@lheapop.gsfc.nasa.gov
*
*/
#ifndef PHOTON_INFO_H
#define PHOTON_INFO_H

class PhotonInfo
{
public:
    void setTime(double time)       { m_time=time; }
    void setEnergy(double energy)   { m_energy=energy; }
    void setTheta(float theta)      { m_theta=theta; }
    void setPhi(float phi)          { m_phi=phi; }
    void setSignal(bool signal)     { m_signal=signal; }
    
    inline double time() const     { return m_time; }
    inline double energy() const   { return m_energy; }
    inline float theta() const     { return m_theta; }
    inline float phi() const       { return m_phi; }
    inline bool signal() const     { return m_signal; }
    
private:
    double m_time;
    double m_energy;
    float  m_theta;
    float  m_phi;
    bool   m_signal;
};

#endif
