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
ISpectrumFactory & GaussianSourceFactory();
ISpectrumFactory & GRBmanagerFactory();
ISpectrumFactory & GRBobsFactory();
ISpectrumFactory & IsotropicFactory();
ISpectrumFactory & MapSourceFactory();
ISpectrumFactory & MapCubeFactory();
ISpectrumFactory & PeriodicSourceFactory();
ISpectrumFactory & PulsarFactory();
ISpectrumFactory & PulsarSpectrumFactory();
ISpectrumFactory & SimpleTransientFactory();
ISpectrumFactory & SpectralTransientFactory();
ISpectrumFactory & TransientTemplateFactory();

SpectrumFactoryLoader::SpectrumFactoryLoader()
{
   std::vector<ISpectrumFactory *> factories;

   factories.push_back(&GRBmanagerFactory());
   factories.push_back(&GRBobsFactory());
   factories.push_back(&GaussianSourceFactory());
   factories.push_back(&IsotropicFactory());
   factories.push_back(&MapSourceFactory());
   factories.push_back(&MapCubeFactory());
   factories.push_back(&PeriodicSourceFactory());
   factories.push_back(&PulsarFactory());
   factories.push_back(&SimpleTransientFactory());
   factories.push_back(&SpectralTransientFactory());
   factories.push_back(&TransientTemplateFactory());
   factories.push_back(&PulsarSpectrumFactory());

//    std::vector<ISpectrumFactory *>::const_iterator factory;
//    for (factory = factories.begin(); factory != factories.end(); ++factory){
//       m_names.push_back(factory->name());
//    }
}
