/*!
* \class BackgroundMixer
*
* \brief This class provides interface needed to create the background mixed events data.
* \author Jay Norris        jnorris@lheapop.gsfc.nasa.gov
* \author Sandhia Bansal    sandhiab@lheapop.gsfc.nasa.gov
*
* This class implements methods to create background mixed data from specified events and
* background data.  Alternatively, it provides methods to read already mixed data.
*
*/


#ifndef BACKGROUND_MIXER_H
#define BACKGROUND_MIXER_H

#include <vector>
#include <string>
#include <fstream>


#include "PhotonInfo.h"


class BackgroundMixer
{
public:
    /*!
     * \brief Constructors.
     */
    BackgroundMixer(const std::string &grbFile, const std::string &backgroundFile, 
        const std::string &mixedDir, const double grbOffsetTime);
    BackgroundMixer(const std::string &mixedFile);
    

    
    // No memory management function required
    // So - need for destructor, copy constructor and assignment operator to be defined
   
    

    // Accessors
    std::string grbFile() const                        { return m_grbFile; }
    std::string backgroundFile() const                 { return m_backgroundFile; }
    std::string mixedDir() const                       { return m_mixedDir; }
    std::string mixedFile() const                      { return m_mixedFile; }
    double grbOffsetTime() const                       { return m_grbOffsetTime; }
    const std::vector<PhotonInfo> &photonData() const  { return m_photonData; }
    long ngrb() const                                  { return m_nGrb; }
    long nbck() const                                  { return m_nBck; }
    
    
private:
    /*!
     * \brief Create background mixed events data.
     */
    void mix();

    /*!
     * \brief Read background data from specfied file.
     */
    void readBackground();

    /*!
     * \brief Read events data from specified file.
     */
    void readGRB();

    /*!
     * \brief Read background mixed data.
     */
    void readMixedData();

    /*!
     * \brief Write background mixed data to a file for use by other programs */
    void writeMixedData();
    

    
    // data members
    std::string              m_grbFile;              //! Name of file containing events data.
    std::string              m_backgroundFile;       //! Name of file that lists background data.
    std::string              m_mixedDir;             //! Directory where mixed file does or will reside.
    std::string              m_mixedFile;            //! File that contains or will contain background mixed data.
    double					 m_grbOffsetTime;        //! Time needed to offset events data times.
    std::vector<PhotonInfo>  m_photonData;           //! Photon data list ((time,energy) pairs).
    long                     m_nGrb;                 //! Number of photons in events data.
    long                     m_nBck;                 //! Number of photons read from background data file.
};



/*!
 * \brief This class implements method to sort data in time order.
 */
class timeCmp
{
public:
    bool operator()(PhotonInfo &data1, PhotonInfo &data2)
    {
        return data1.time() < data2.time();    
    }
};



#endif // BACKGROUND_MIXER_H
