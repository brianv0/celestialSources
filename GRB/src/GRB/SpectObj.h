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
 *
 * \author Nicola Omodei       nicola.omodei@pi.infn.it 
 * \author Johann Cohen-Tanugi johann.cohen@pi.infn.it
 *
 */
class SpectObj
{
 public:
  /// empty default constructor necessary to make a SpectObj instance 
  /// a data member of other classes
  SpectObj(){;}

  /*! Constructor using values for defining the set of energies 
   * in logarithmic scale
   * @param enmin minimal energy
   * @param enmax maximal energy
   * @param enstep number of steps from enmin to enmax 
   */
   SpectObj(double,double,double);


   /*! Constructor using an already existing map <energy,value>
    * @param begin const_iterator pointing to the desired start
    * @param end const_iterator pointing to the desired end
    */
  SpectObj(std::map<double,double>::const_iterator start,
	   std::map<double,double>::const_iterator end);

  ~SpectObj(){;}
  
  /*! compute \f$\sum E_ip(E_i)B_i\f$, where \f$ p(E_i)\f$ is 
   * the spectrum value at energy \f$ i\f$ and \f$B_i\f$ is the bin size of 
   * the bin starting at \f$ i\f$
   * @param emin minimal energy \f$E_i\f$
   * @param emin maximal energy \f$E_i\f$
   */
  double integrated_E_Flux(double emin, double emax);
  
  /*! compute \f$\sum p(E_i)B_i\f$, where \f$p(E_i)\f$ is 
   * the spectrum value at energy \f$ i\f$ and \f$B_i\f$ is the bin size of 
   * the bin starting at \f$ i\f$
   * @param emin minimal energy \f$E_i\f$
   * @param emin maximal energy \f$E_i\f$
   */
  double integrated_Flux(double emin, double emax);
  
  /*! returns the size value of the bin starting at it
   * \param it const iterator pointing to the min energy of the bin 
   */
   double binSize(std::map<double,double>::const_iterator);
  
  /// returns the size of the map, aka the number of values points on the spectrum
  inline double size()  { return m_spectrum.size(); }

  inline void   clear() { m_spectrum.clear();       }

  inline double getBinContent(double index) 
    {
      std::map<double,double>::const_iterator it = m_spectrum.begin();
      std::advance(it,index);
      return it->second;
    }

  /*! returns the i-th energy value
   * \param i index in the map of the desired energy 
   */
  inline double getBinValue(double index) 
    {
      std::map<double,double>::const_iterator it = m_spectrum.begin();
      std::advance(it,index);
      return it->first;
    }
  void clearSpectrumVector();
  //////////////////////////////////////////////////
  /// adds a spectrum object
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

  SpectObj operator/(const double value);

  inline SpectObj operator/=(const double value)
  {
    return (*this)/value;
  }
  //////////////////////////////////////////////////

  /*! returns vector of bins, aka of values E_i
   * @param value scale factor
   */
  std::vector<double> getEnergyVector(const double value = 1.);

  /*! returns vector of bins, aka of values \f$E_{i+1} - E_i\f$
   * @param value scale factor
   */
  std::vector<double> getBinVector(const double value = 1.);

  std::vector<double> getSpectrumVector(const double value = 1.);
  
  /*! creates a new SpectObj with energy range [start,end]
   * \param start const iterator positioned at the min of the energy range
   * \param end  const iterator positioned at the max of the energy range
   */
  SpectObj extractSub(std::map<double,double>::const_iterator /*start*/,
                      std::map<double,double>::const_iterator /*end*/);

  /*! Change one spectral value
   * @param energy position of value to be changed
   * @param new value
   */
  SpectObj SetSpectrum(double energy, double spectrum);

  /*! Load this SpectObj instance with energy and spectrum vectors
   * @param energy vector of energy values
   * @param spectrum vector of spectral values
   */
  SpectObj SetSpectrum(std::vector<double> energy, std::vector<double> spectrum);

  SpectObj AddSpectrum(SpectObj);

 private:
  /// internal representation 
  std::map<double,double> m_spectrum;
  
};
#endif
