/**
 *
 */

#include "genericSources/FILESpectrumMap.h"
#include "flux/SpectrumFactory.h"
#include "facilities/Util.h"
#include <iostream>

ISpectrumFactory &FILESpectrumMapFactory() {
  static SpectrumFactory<FILESpectrumMap> factory;
  return factory;
}


FILESpectrumMap::FILESpectrumMap(const std::string& params)
  : MapSource(params)
{
  facilities::Util::keyValueTokenize(params,",",m_parmap);

  std::string file_name = m_parmap["fileName"];
  facilities::Util::expandEnvVar(&file_name);
  std::string islog=m_parmap["log"];

  double e_min = std::atof(m_parmap["emin"].c_str());
  double e_max = std::atof(m_parmap["emax"].c_str());

  std::string filespectrum_params=file_name;
  if(islog.find("yes")!=std::string::npos){
    filespectrum_params+=",log";
  }
  m_filespectrum = new FILESpectrum(filespectrum_params);

  if(m_flux==0.)
    m_flux = m_filespectrum->flux(0);
}


