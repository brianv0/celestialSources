/** 
* @file SpectrumFactoryLoader.cxx
* @brief Load the external spectrum factory objects
*
*  $Header$
*/

#include <iostream>
#include <vector>

#include "celestialSources/SpectrumFactoryLoader.h"
#ifndef BUILD_WITHOUT_ROOT
#include "celestialSources/TRandom4.h"
#endif

#include "flux/ISpectrumFactory.h"

// declare the external factories.
ISpectrumFactory & FitsTransientFactory();
ISpectrumFactory & GaussianSourceFactory();
ISpectrumFactory & GRBmanagerFactory();
ISpectrumFactory & GRBobsmanagerFactory();
ISpectrumFactory & GRBtemplateManagerFactory();
ISpectrumFactory & IsotropicFactory();
ISpectrumFactory & MapSourceFactory();
ISpectrumFactory & MapCubeFactory();
ISpectrumFactory & PeriodicSourceFactory();
ISpectrumFactory & PulsarFactory();
ISpectrumFactory & PulsarSpectrumFactory();
ISpectrumFactory & SourcePopulationFactory();
ISpectrumFactory & SimpleTransientFactory();
ISpectrumFactory & SpectralTransientFactory();
ISpectrumFactory & TransientTemplateFactory();
ISpectrumFactory & TF1SpectrumFactory();
ISpectrumFactory & TF1MapFactory();
ISpectrumFactory & FileSpectrumFactory();
ISpectrumFactory & FileSpectrumMapFactory();
ISpectrumFactory & microQuasarFactory();
ISpectrumFactory & RadialSourceFactory();

SpectrumFactoryLoader::SpectrumFactoryLoader() {
// Replace ROOT's global TRandom instance with our local version that
// uses the CLHEP engines underneath.
#ifndef BUILD_WITHOUT_ROOT
   gRandom = new TRandom4();
#endif

   load(FitsTransientFactory());
#ifndef BUILD_WITHOUT_ROOT
   load(GRBmanagerFactory());
   load(GRBobsmanagerFactory());
   load(GRBtemplateManagerFactory());
#endif
   load(GaussianSourceFactory());
   load(IsotropicFactory());
   load(MapSourceFactory());
   load(MapCubeFactory());
   load(PeriodicSourceFactory());
   load(PulsarFactory());
   load(SimpleTransientFactory());
   load(SourcePopulationFactory());
   load(SpectralTransientFactory());
   load(TransientTemplateFactory());
#ifndef BUILD_WITHOUT_ROOT
   load(PulsarSpectrumFactory());
   load(TF1SpectrumFactory());
   load(TF1MapFactory());
#endif
   load(FileSpectrumFactory());
   load(FileSpectrumMapFactory());
   load(microQuasarFactory());
   load(RadialSourceFactory());
}

void SpectrumFactoryLoader::load(ISpectrumFactory& factory) {
       m_names.push_back( factory.name() );
}
