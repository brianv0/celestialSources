/**
  * @file TF1Map.cxx
  * @author Johann Cohen-Tanugi
  * $Header$
  */

#include "genericSources/TF1Map.h"
#include "flux/SpectrumFactory.h"
#include "facilities/Util.h"
#include <cstdlib>
#include <iostream>

ISpectrumFactory &TF1MapFactory() {
  static SpectrumFactory<TF1Map> factory;
  return factory;
}


TF1Map::TF1Map(const std::string& params)
  : MapSource(params)
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
