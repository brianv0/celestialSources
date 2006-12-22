/**
 * @file test.cxx
 * @brief Code to check SpectrumFactoryLoader class
 * $Header$
 */

#include <iostream>
#include <string>
#include <vector>

#include "celestialSources/SpectrumFactoryLoader.h"

int main() {
   
   SpectrumFactoryLoader loader;

   const std::vector<std::string> & names = loader.names();
   std::cout << "Loaded the following factories: \n";
   for (unsigned int i = 0; i < names.size(); i++) {
      std::cout << names.at(i) << "\n";
   }

}
