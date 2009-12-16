/**
 * @file FileSpectrumMap.cxx
 * @author Johann Cohen-Tanugi
 * $Header$
 */

#include "genericSources/FileSpectrumMap.h"
#include "flux/SpectrumFactory.h"
#include "facilities/Util.h"
#include <cstdlib>
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

  //MapSource will not set e_min and e_max properly if
  //there is no gamma defined. The powerlaw index being irrelevant
  //here, I just read the parameters here...
  m_emin = std::atof(m_parmap["emin"].c_str());
  m_emax = std::atof(m_parmap["emax"].c_str());
}


