/** 
* @file SpectrumFactoryLoader.h
* @brief Load the external spectrum factory objects
*
*  $Header$
*/

#ifndef SpectrumFactoryLoader_h
#define SpectrumFactoryLoader_h

#include <vector>
#include <string>

class SpectrumFactoryLoader {
public:
    /// ctor does the work
    SpectrumFactoryLoader();
    /// access to a list of the names that were loaded
    std::vector<std::string> names()const{return m_names;}
private:
    std::vector<std::string> m_names;
};

#endif
