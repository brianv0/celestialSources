#include <iterator>
#include <map>
#include <vector>
#include <math.h>
#ifndef SpectObj_H
#define SpectObj_H
/*!
 * \class SpectObj
 *
 * \brief This class describe the Spectrum Object.
 * 
 * A Spectrum Object is basically a map containing in the first entry 
 * the energy and in the second entry the value of the spectrum.
 * The spectrum is in \f$ph/s/MeV\f$
 *
 * \author Nicola Omodei       nicola.omodei@pi.infn.it 
 * \author Johann Cohen-Tanugi johann.cohen@pi.infn.it
 *
 */
class SpectObj
{
 public:
  SpectObj(){;}

  /*!Constructor using values for defining the set of energies 
   * in logarithmic scale
   *
   * \param enmin minimal energy
   * \param enmax maximal energy
   * \param enstep number of steps from enmin to enmax 
   */
  SpectObj(double enmin, double enmax ,double enstep);
  
  
  /*!Constructor using an already existing map <energy,value>
   * \param start const_iterator pointing to the desired start
   * \param end const_iterator pointing to the desired end
   */
  SpectObj(std::map<double,double>::const_iterator start,
	   std::map<double,double>::const_iterator end);

  ~SpectObj(){;}
  
  /*!Computes \f$\sum E_i p(E_i) B_i\f$, where \f$p(E_i)\f$ is 
   * the spectrum value at energy \f$E_i\f$ and \f$B_i\f$ is the bin size of 
   * the bin starting at \f$i\f$
   * \param enmin minimal energy \f$E_i\f$
   * \param enmax maximal energy \f$E_i\f$
   * \retval eV/s
   */
  double integrated_E_Flux(double enmin, double enmax);
  
  /*! computes \f$\sum p(E_i)B_i\f$, where \f$p(E_i)\f$ is 
   * the spectrum value at energy \f$E_i\f$ and \f$B_i\f$ is the bin size of 
   * the bin starting at \f$i\f$
   * \param enmin minimal energy \f$E_i\f$
   * \param enmax maximal energy \f$E_i\f$
   * \retval ph/s
   */
  double integrated_Flux(double enmin, double enmax);
  
  /*! returns the size value of the bin starting at it
   * \param it const iterator pointing to the min energy of the bin 
   * \retval eV
   */
  double binSize(std::map<double,double>::const_iterator it);
  
  /*! returns the size of the map, aka the number of values points on the spectrum
   */
  inline double size()  { return m_spectrum.size(); }
  /*! It empties the values of the spectrum, not the energy vector!
   */
  inline void   clear() { m_spectrum.clear();       }
   
  /*! returns the spectrum in the i-th energy value
   * \param index is the position in the map of the desired energy 
   * \retval eV/MeV/s
   */
  inline double getBinContent(double index) 
    {
      std::map<double,double>::const_iterator it = m_spectrum.begin();
      std::advance(it,index);
      return it->second;
    }
  
  /*! returns the i-th energy value
   * \param index is the position in the map of the desired energy 
   * \retval eV
   */
  inline double getBinValue(double index) 
    {
      std::map<double,double>::const_iterator it = m_spectrum.begin();
      std::advance(it,index);
      return it->first;
    }
  
  //////////////////////////////////////////////////
  //! adds a spectrum object
  SpectObj operator+(const SpectObj inputObj);
  
  inline SpectObj operator+=(const SpectObj inputObj)
  {
    return (*this)+inputObj;
  }

  /// scales the spectrum values
  SpectObj operator*(const double value);

  inline SpectObj operator*=(const double value)
  {
    return (*this)*value;
  }
  //! divide a spectrum object for a value
  SpectObj operator/(const double value);
  
  inline SpectObj operator/=(const double value)
  {
    return (*this)/value;
  }
  //////////////////////////////////////////////////

  /*! returns vector of bins, aka of values \f$E_i\f$
   * \param value scale factor
   */
  std::vector<double> getEnergyVector(const double value = 1.);
  
  /*! returns vector of bins, aka of values \f$E_{i+1} - E_i\f$
   * \param value scale factor
   */
  std::vector<double> getBinVector(const double value = 1.);
  /*! returns a vector containing the values of the flux \f$p(E_i)\f$
   * \param value scale factor
   */
  std::vector<double> getSpectrumVector(const double value = 1.);
  
  /*! creates a new SpectObj with energy range [start,end]
   * \param start const iterator positioned at the min of the energy range
   * \param end  const iterator positioned at the max of the energy range
   */
  SpectObj extractSub(std::map<double,double>::const_iterator start,
		      std::map<double,double>::const_iterator end);
  
  /*! Changes one spectral value
   * \param energy position of value to be changed
   * \param spectrum new value
   */
  SpectObj SetSpectrum(double energy, double spectrum);
  
  /*! Loads this SpectObj instance with energy and spectrum vectors
   * \param energy vector of energy values
   * \param spectrum vector of spectral values
   */
  SpectObj SetSpectrum(std::vector<double> energy, std::vector<double> spectrum);
  
  
 private:
  /// internal representation 
  std::map<double,double> m_spectrum;
  
};
#endif
