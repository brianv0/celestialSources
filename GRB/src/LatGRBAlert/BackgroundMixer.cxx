// FILE: BackgroundMixer.cxx
// This class reads GRB and background photons lists and mixes them together.


#include "BackgroundMixer.h"
#include <iostream>
#include <algorithm>              // for sort
#include <iomanip>
#ifndef WIN32
#include <stdio.h>
#endif
#include "facilities/Util.h"




// Constructor  BackgroundMixer()
// Read GRB and background photon lists, mix and write the mixed list to an output file.
BackgroundMixer::BackgroundMixer(const std::string &grbFile, 
                                 const std::string &backgroundFile, 
                                 const std::string &mixedFile, 
                                 const double grbOffsetTime)
               :m_grbFile(grbFile),
                m_backgroundFile(backgroundFile),
                m_mixedFile(mixedFile),
                m_grbOffsetTime(grbOffsetTime),
                m_photonData()
{
    readGRB();
    readBackground();
    mix();
    
    writeMixedData();
}



// Constructor  BackgroundMixer()
// Reads an already mixed photon list.
BackgroundMixer::BackgroundMixer(const std::string &mixedFile)
               :m_mixedFile(mixedFile),
                m_photonData()
{
    readMixedData();
}



// readGRB()
// Reads a file containing GRB photon list.
void BackgroundMixer::readGRB()
{
    std::cout << "read " << m_grbFile << std::endl;
    try 
    {
        facilities::Util::expandEnvVar(&m_grbFile);
        std::ifstream is(m_grbFile.c_str());
        
        if (is)
        {
            std::string str;
            is >> str >> m_nGrb;
            
            m_photonData.resize(m_nGrb);
            
            double t, e, th, p;
            for (long i=0; i<m_nGrb; ++i)
            {
                is >> t >> e >> th >> p;
                m_photonData[i].setTime(t+m_grbOffsetTime);
                m_photonData[i].setEnergy(e);
                m_photonData[i].setTheta(th);
                m_photonData[i].setPhi(p);
                m_photonData[i].setSignal(1);
            }
        }
        
        else
            std::cout << "Error while opening " << m_grbFile << std::endl;
    }
    
    catch (...)
    {
        std::cout << "Error while reading " << m_grbFile << std::endl;
    }
}



// readBackground()
// Reads a file containing background photon list.
// for now assume that background file has same structure as the GRB file
void BackgroundMixer::readBackground()
{
    std::cout << "read " << m_backgroundFile << std::endl;
    try 
    {
        facilities::Util::expandEnvVar(&m_backgroundFile);
        std::ifstream is(m_backgroundFile.c_str());
        
        if (is)
        {
            std::string str;
            is >> str >> m_nBck;
            
            std::vector<PhotonInfo>::size_type sz = m_nGrb + m_nBck;
            m_photonData.resize(sz);
            
            double t, e, th, p;
            for (long i=m_nGrb; i<sz; ++i)
            {
                is >> t >> e >> th >> p;
                m_photonData[i].setTime(t);
                m_photonData[i].setEnergy(e);
                m_photonData[i].setTheta(th);
                m_photonData[i].setPhi(p);
                m_photonData[i].setSignal(0);
            }
        }
        
        else
            std::cout << "Error while opening " << m_backgroundFile << std::endl;
    }
    
    catch (...)
    {
        std::cout << "Error while reading " << m_backgroundFile << std::endl;
    }
}



// mix()
// Reads a file containing mixed photon list.
void BackgroundMixer::mix()
{
    std::sort(m_photonData.begin(), m_photonData.end(), timeCmp());
}



// writeMixedData()
// Writes mixed photon list to an output file
void BackgroundMixer::writeMixedData()
{
    facilities::Util::expandEnvVar(&m_mixedFile);
    FILE *fptr = fopen(m_mixedFile.c_str(), "w");
    
    if (m_photonData.size() > 0)
    {
        fprintf(fptr, "nGRB: %d   nBCK: %d \n", m_nGrb, m_nBck);
        for (std::vector<PhotonInfo>::iterator it=m_photonData.begin(); it!=m_photonData.end(); ++it)
            fprintf(fptr,"%12.5lf\t%12.5lf\t%12.5f\t%12.5f\t%d\n", 
            it->time(), it->energy(), it->theta(), it->phi(), it->signal());
    }
    
    fclose(fptr);
}



// readMixedData()
// Reads mixed photon list from a file
void BackgroundMixer::readMixedData()
{
    std::cout << "read " << m_mixedFile << std::endl;
    try 
    {
        facilities::Util::expandEnvVar(&m_mixedFile);
        std::ifstream is(m_mixedFile.c_str());
        
        if (is)
        {
            std::string  str;
            
            is >> str >> m_nGrb >> str >> m_nBck;
            
            long sz = m_nGrb + m_nBck;
            
            m_photonData.resize(sz);
            
            double t, e, th, ph;
            bool   s;
            for (long i=0; i<sz; ++i)
            {
                is >> t >> e >> th >> ph >> s;
                m_photonData[i].setTime(t);
                m_photonData[i].setEnergy(e);
                m_photonData[i].setTheta(th);
                m_photonData[i].setPhi(ph);
                m_photonData[i].setSignal(s);
            }
        }
        
        else
            std::cout << "Error while opening " << m_mixedFile << std::endl;
    }
    
    catch (...)
    {
        std::cout << "Error while reading " << m_mixedFile << std::endl;
    }
}


