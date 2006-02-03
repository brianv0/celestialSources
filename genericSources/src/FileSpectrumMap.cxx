/**
 *
 */

#include "genericSources/FileSpectrumMap.h"
#include "flux/SpectrumFactory.h"
#include "facilities/Util.h"
#include <iostream>

ISpectrumFactory &FileSpectrumMapFactory() {
  static SpectrumFactory<FileSpectrumMap> factory;
  return factory;
}


FileSpectrumMap::FileSpectrumMap(const std::string& params)
  : MapSource(params)
{
  m_filespectrum = new FileSpectrum(params);
  //Let FileSpectrum decide which flux to use
  m_flux = m_filespectrum->flux();

  facilities::Util::keyValueTokenize(params,",",m_parmap);

  double e_min = std::atof(m_parmap["emin"].c_str());
  double e_max = std::atof(m_parmap["emax"].c_str());

}


