/**  @file test_main.cxx
@brief microQuasar test program

$Header$
*/

#include "flux/SpectrumFactoryTable.h"
#include "flux/FluxMgr.h"

#include <iostream>
#include <algorithm>


namespace {
    int default_count =100 ;
    std::string default_source("Galactic_LSI61_303");
    std::string xml_spec("$(MICROQUASARROOT)/xml/mq.xml");

}
#define DECLARE_SPECTRUM(x)   extern const ISpectrumFactory& x##Factory; x##Factory.addRef();


void listSources(const std::list<std::string>& source_list ) {
    std::cout << "List of available sources:" << std::endl;
    for( std::list<std::string>::const_iterator it = source_list.begin(); 
        it != source_list.end(); ++it) { 
            std::cout << '\t'<< *it << std::endl;
        }

}
void listSpectra() {
    std::cout << "List of loaded Spectrum objects: " << std::endl;
    std::list<std::string> spectra(SpectrumFactoryTable::instance()->spectrumList());
    for( std::list<std::string>::const_iterator it = spectra.begin(); 
        it != spectra.end(); ++it) { 
            std::cout << '\t'<< *it << std::endl;
        }
}




int main(int argn, char * argc[]) {
    using std::cout;
    using std::endl;
    int rc(0);

    try {
        // load our factory by name
        DECLARE_SPECTRUM( microQuasar);

        int count = default_count;
        std::string source_name(default_source);

        // set up the flux manager, adding in our xml description
        std::vector<std::string> fileList;
        fileList.push_back("$(FLUX_XML)/source_library.xml");
        fileList.push_back(xml_spec);
        FluxMgr fm(fileList);


        if ( argn >1 ) source_name = argc[1];
        if( source_name =="list") { 
            listSources(fm.sourceList());
            listSpectra(); return 0; }
        if ( argn >2 ) count = ::atoi(argc[2]);

        cout << "------------------------------------------------------" << endl;
        cout << " MQ test: " << endl;
        cout << ( ( argn ==1)?  " No command line args, using default flux \""
            :  " Selected source name \"");
        cout  << source_name <<"\"" << endl;


        std::list<std::string> source_list(fm.sourceList());
        if( std::find(source_list.begin(), source_list.end(), source_name)==source_list.end() ) {

            std::cout << "Source \"" << source_name << "\" not found in the list of sources!" << std::endl;
            listSources(source_list);
            return -1;
        }

        fm.test(std::cout, source_name, count); 

    }catch(const std::exception & x){
        std::cerr << "Caught exception: " << x.what() << std::endl;
        rc =1 ;
    }

    return rc;    
}

