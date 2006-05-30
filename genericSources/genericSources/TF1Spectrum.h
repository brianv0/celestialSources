/**
 * @file TF1Spectrum.h
 * @brief A source class that uses ROOT TF1 object to define a spectrum symbolically
 * id est with a formula. The class is instantiated based on the information parsed from
 * the params string, passed to the object vi the constructor. Here is an example of
 * such params string :
 * params="flux=17.,tf1name=TF1Spectrum_TEST,formula=-0.0001*(100.-x)*(1100.-x),emin=100.,emax=1100"
 * tf1name is needed to let the user make sure that the ROOT object has a unique name identifier in ROOT internal memory.
 * tf1precision defines the granularity with which TF1 will approximate the formula.
 * @author Johann Cohen-Tanugi
 *
 * $Header$
 */
 
#ifndef TF1SPECTRUM_H
#define TF1SPECTRUM_H

#include "TF1.h"
#include "flux/Spectrum.h"
#include<map>

class TF1Spectrum : public Spectrum
{
 public:
  TF1Spectrum(const std::string& /*params*/);
  ~TF1Spectrum() {;}
  
  ///The TF1 object manages the random draw, so that the float argument, normally a
  /// random draw from a [0,1[ distribution, is unused.
  float operator()(float )const 
    {
      return p_tf1.GetRandom();
    }

  std::string title() const 
    {
      return "TF1Spectrum";
    }

  const char * particleName() const
    {
      return m_particle_name.c_str();
    }

  
  ///Return an energy sampled from the TF1 distribution. 
  double energy(double time)
    {
      return (*this)(time);
    }

 private:
  ///pointer to the TF1 object that reads the symbolic formula and does all the job
  mutable TF1  p_tf1;
  std::map<std::string,std::string> m_parmap;
};

#endif
