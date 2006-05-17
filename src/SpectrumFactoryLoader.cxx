/** 
* @file SpectrumFactoryLoader.cxx
* @brief Load the external spectrum factory objects
*
*  $Header$
*/

#include <iostream>
#include <vector>

#include "celestialSources/SpectrumFactoryLoader.h"
#include "celestialSources/TRandom4.h"

#include "flux/ISpectrumFactory.h"

// declare the external factories.
ISpectrumFactory & FitsTransientFactory();
ISpectrumFactory & GaussianSourceFactory();
ISpectrumFactory & GRBmanagerFactory();
ISpectrumFactory & GRBobsmanagerFactory();
ISpectrumFactory & GRBobsFactory();
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

SpectrumFactoryLoader::SpectrumFactoryLoader() {
// Replace ROOT's global TRandom instance with our local version that
// uses the CLHEP engines underneath.
   gRandom = new TRandom4();

   load(FitsTransientFactory());
   load(GRBmanagerFactory());
   load(GRBobsmanagerFactory());
   load(GRBobsFactory());
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
   load(PulsarSpectrumFactory());
   load(TF1SpectrumFactory());
   load(TF1MapFactory());
   load(FileSpectrumFactory());
   load(FileSpectrumMapFactory());
}

void SpectrumFactoryLoader::load(ISpectrumFactory& factory) {
       m_names.push_back( factory.name() );
}
