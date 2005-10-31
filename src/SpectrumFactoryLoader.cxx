/** 
* @file SpectrumFactoryLoader.cxx
* @brief Load the external spectrum factory objects
*
*  $Header$
*/

#include <vector>

#include "celestialSources/SpectrumFactoryLoader.h"
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

SpectrumFactoryLoader::SpectrumFactoryLoader()
{
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
}
void SpectrumFactoryLoader::load(ISpectrumFactory& factory)
{
       m_names.push_back( factory.name() );
}

