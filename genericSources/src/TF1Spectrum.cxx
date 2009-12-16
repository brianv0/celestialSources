/**
 * @file TF1Spectrum.cxx
 * @author Johann Cohen-Tanugi
 * $Header$
 */

#include "genericSources/TF1Spectrum.h"
#include "flux/SpectrumFactory.h"
#include "facilities/Util.h"
#include <cstdlib>
#include <iostream>

ISpectrumFactory &TF1SpectrumFactory() {
  static SpectrumFactory<TF1Spectrum> factory;
  return factory;
}


TF1Spectrum::TF1Spectrum(const std::string& params)
  : Spectrum()
{
  facilities::Util::keyValueTokenize(params," \t,\n",m_parmap);
  std::string internal_name = m_parmap["tf1name"].c_str();
  int grid_bins = 0;
  grid_bins = std::atoi(m_parmap["tf1precision"].c_str());

  m_flux = std::atof(m_parmap["flux"].c_str());
  double e_min = std::atof(m_parmap["emin"].c_str());
  double e_max = std::atof(m_parmap["emax"].c_str());

  p_tf1 = TF1(internal_name.c_str(),m_parmap["formula"].c_str(), e_min, e_max);

  if(grid_bins!=0){
    p_tf1.SetNpx(grid_bins);
  }

  if(m_flux==0.)
    m_flux = p_tf1.Integral(e_min,e_max);

}


