// $Header$

#include "../../GRBmaker/GRBmaker.h"
#include "src/FluxMgr.h"
//#include "flux/EventSource.h"
//#include "flux/ISpectrumFactory.h"
//#include "../SpectrumFactoryTable.h"
#include <iostream>
#include <string>
//#include <algorithm>
//#include "flux/CHIMESpectrum.h"
//#include "flux/Orbit.h"


//static int default_count = 10;
//Testing
//static const char * default_source="default";
//Default
//static const char * default_source="CrElectron";

//void help() {
//    std::cout << 
//        "   Simple test program to create a particle source, then run it.\n"
//        "   Command line args are \n"
//        "      <source name> <count>\n"
//        "   with defaults \"" 
//        <<  default_source << "\"," << default_count
//        << "\n  Also, 'help' for this help, 'list' for a list of sources and Spectrum objects "
//        << std::endl;
//}
//void listSources(const std::list<std::string>& source_list ) {
//    std::cout << "List of available sources:" << std::endl;
//    for( std::list<std::string>::const_iterator it = source_list.begin(); 
//    it != source_list.end(); ++it) { 
//        std::cout << '\t'<< *it << std::endl;
//    }
//    
//}
//void listSpectra() {
//    std::cout << "List of loaded Spectrum objects: " << std::endl;
//    std::list<std::string> spectra(SpectrumFactoryTable::instance()->spectrumList());
//    for( std::list<std::string>::const_iterator it = spectra.begin(); 
//    it != spectra.end(); ++it) { 
//        std::cout << '\t'<< *it << std::endl;
//   }
//}


//#define DLL_DECL_SPECTRUM(x)   extern const ISpectrumFactory& x##Factory; x##Factory.addRef();

//void flux_load() {

// these are the spectra that we want to make available
//    DLL_DECL_SPECTRUM( CHIMESpectrum);
//    DLL_DECL_SPECTRUM( AlbedoPSpectrum);
//    DLL_DECL_SPECTRUM( HeSpectrum);
//    DLL_DECL_SPECTRUM( GalElSpectrum);
//    DLL_DECL_SPECTRUM( CrElectron);
//    DLL_DECL_SPECTRUM( CrProton);
//  DLL_DECL_SPECTRUM( CREMESpectrum);
//}

//int main(int argn, char * argc[]) {
int main()
{
    GRBmaker *grbMaker(new GRBmaker);
    grbMaker->create(36.804426, 18, 9.843296, 0.199746, 0.564782, 2.07933, 2.40867, 0.000351, 1);

    delete grbMaker;
    
    
    return 0;
}

void WARNING (const char * text ){  std::cerr << "WARNING: " << text << '\n';}
void FATAL(const char* s){std::cerr << "\nERROR: "<< s;}
