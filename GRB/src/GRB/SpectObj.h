#ifndef SpectObj_H
#define SpectObj_H
#include <iterator>
#include <map>
#include <vector>
#include <cmath>
#include "GRBConstants.h"
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
  double integrated_E_Flux(double enmin, double enmax=cst::enmax) const;
  
  /*! computes \f$\sum p(E_i)B_i\f$, where \f$p(E_i)\f$ is 
   * the spectrum value at energy \f$E_i\f$ and \f$B_i\f$ is the bin size of 
   * the bin starting at \f$i\f$
   * \param enmin minimal energy \f$E_i\f$
   * \param enmax maximal energy \f$E_i\f$
   * \retval ph/s
   */
  double integrated_Flux(double enmin, double enmax=cst::enmax) const;
  
  /*! returns the size value of the bin starting at it
   * \param it const iterator pointing to the min energy of the bin 
   * \retval eV
   */
  double binSize(std::map<double,double>::const_iterator it) const;
  
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
  inline double getSpectrum(int index) const
    {
      std::map<double,double>::const_iterator it = m_spectrum.begin();
      std::advance(it,index);
      return it->second;
    }
  
  /*! returns the i-th energy value
   * \param index is the position in the map of the desired energy 
   * \retval eV
   */
  inline double getEnergy(int index) const
    {
      std::map<double,double>::const_iterator it = m_spectrum.begin();
      std::advance(it,index);
      return it->first;
    }
  
  inline double getBin(int index) const 
    {
      double e1,e2;
      std::map<double,double>::const_iterator it = m_spectrum.begin();
      std::advance(it,index);
      e1 = it->first;
      std::advance(it,1);
      e2 = it->first;
      return e2-e1;
    }
  
  //////////////////////////////////////////////////
  //! adds a spectrum object
  const SpectObj operator+(const SpectObj inputObj);
  
  inline SpectObj operator+=(SpectObj inputObj)
  {
    return (*this)+inputObj;
  }

  /// scales the spectrum values
  const SpectObj operator*(const double value);
  inline SpectObj operator*=(const double value)
  {
    return (*this)*value;
  }
  //! divide a spectrum object for a value
  const SpectObj operator/(const double value);
  
  inline SpectObj operator/=(const double value)
  {
    return (*this)/value;
  }
  
  
  SpectObj em2obs(double gamma, double angle);
  SpectObj obs2em(double gamma, double angle);
  inline double SpectObj::getEnergyTransformation(double gamma, double angle)
    {
      double beta = sqrt(1.-1./pow(gamma,2.));
      return gamma*(1+beta*cos(angle));
    }
  //////////////////////////////////////////////////

  /*! returns vector of bins, aka of values \f$E_i\f$
   * \param value scale factor
   */
  std::vector<double> getEnergyVector(const double value = 1.) const;
  
  /*! returns vector of bins, aka of values \f$E_{i+1} - E_i\f$
   * \param value scale factor
   */
  std::vector<double> getBinVector(const double value = 1.) const;
  /*! returns a vector containing the values of the flux \f$p(E_i)\f$
   * \param value scale factor
   */
  std::vector<double> getSpectrumVector(const double value = 1.) const;
  
  /*! creates a new SpectObj with energy range [start,end]
   * \param start const iterator positioned at the min of the energy range
   * \param end  const iterator positioned at the max of the energy range
   */
  SpectObj extractSub(std::map<double,double>::const_iterator start,
		      std::map<double,double>::const_iterator end);
  
  /*! Changes one spectral value
   * \param index position of value to be changed
   * \param spectrum new value
   */
  SpectObj SetSpectrum(int index, double spectrum);
  
  /*! Loads this SpectObj instance with energy and spectrum vectors
   * \param energy vector of energy values
   * \param spectrum vector of spectral values
   */
  SpectObj SetSpectrum(std::vector<double> energy, std::vector<double> spectrum);
  
  /*!
   * \brief returns a photon energy sampled from the current spectrum vector.
   *
   * This is essentially a copycat from the ROOT TH1::GetRandom() method.
   * \param spctrmVec the current spectrum vector \c m_spectrum.
   * \param u uniform random number.
   * \param enmin minimal energy below which no photon energy is drawn. This 
   * is to avoid generation of low energy photons without interest to GLAST.
   */
  
  double DrawPhotonFromSpectrum(float  u = 0.0, 
				double enmin=cst::enmin) const;
  
  
  
 private:
  /// internal representation 
  std::map<double,double> m_spectrum;
  
};
#endif
