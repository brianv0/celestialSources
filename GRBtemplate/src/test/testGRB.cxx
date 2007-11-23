// $Header$

#include "flux/EventSource.h"
#include "flux/SpectrumFactoryTable.h"
#include "flux/FluxMgr.h"
#include "astro/GPS.h"

#include "facilities/commonUtilities.h"

#include <iostream>
#include <fstream>
#include <algorithm>

#include "TTree.h"
#include "TFile.h"

ISpectrumFactory & GRBtemplateManagerFactory();

static int default_count = 10 ;
//Testing
static const char * default_source="GRB_00000";

void help() {
    std::cout << 
        "   Simple test program to create a particle source, then run it.\n"
        "   Command line args are \n"
        "      <source name> <count>\n"
        "   with defaults \"" 
        <<  default_source << "\"," << default_count
        << "\n  Also, 'help' for this help, 'list' for a list of sources and Spectrum objects "
        << std::endl;
}

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


#define DECLARE_SPECTRUM(x)   extern const ISpectrumFactory& x##Factory; x##Factory.addRef();

void flux_load() 
{
  GRBtemplateManagerFactory();
}

void galacticTest(FluxMgr* fm, std::string sourceName)
{
  TTree tree("Events","Events");
  double energy;
  double time=0;
  double l,b,ra,dec,theta,phi;
  tree.Branch("Time",&time,"Time/D");
  tree.Branch("Energy",&energy,"Energy/D");
  tree.Branch("L",&l,"L/D");
  tree.Branch("B",&b,"B/D");
  tree.Branch("Ra",&ra,"Ra/D");
  tree.Branch("Dec",&dec,"Dec/D");
  tree.Branch("Theta",&theta,"Theta/D");
  tree.Branch("Phi",&phi,"Phi/D");

  using astro::GPS;
  EventSource* e = fm->source(sourceName);
  time=fm->time();
  EventSource* f;
  double lavg=0,bavg=0;
  int i=0;
  double interval=0.0;
  while (time<1e8)
    {
      f = e->event(time);
      interval=e->interval(time);
      //here we increment the "elapsed" time and the "orbital" time,
      //just as is done in flux.  NOTE: this is important for the operation 
      //of fluxsource, and is expected.
      time+=interval;
      fm->pass(interval);
      Hep3Vector GlastDir = -(f->launchDir());

      theta = GlastDir.theta();//-TMath::Pi();
      phi   = GlastDir.phi();
      
      Hep3Vector abc(fm->transformToGlast(time,GPS::CELESTIAL).inverse()*(GlastDir));
      astro::SkyDir skyDir(abc,astro::SkyDir::EQUATORIAL);
      
      energy=f->energy();

      l = skyDir.l();
      b = skyDir.b();
      ra = skyDir.ra();
      dec = skyDir.dec();
      if(energy>0.1 && time<1e8)
	{
	  i++;
	  if(i<100) std::cout << "time= "<<time<<" energy= " << energy << " (theta,phi)= " <<theta*180./M_PI<<", "<<phi*180./M_PI<<" interval "<<interval<<std::endl;
	  
	  lavg +=l;
	  bavg +=b;
	  tree.Fill();
	}
    }
  lavg /= i;
  bavg /= i;
  std::cout << "  the average photon location was (l,b) = " << lavg << "," << bavg << std::endl;
  TFile file("GRBTree.root","RECREATE");
  tree.Write();
  file.Close();
  delete e;
}

int main(int argn, char * argc[]) {
  using std::cout;
  using std::endl;
  facilities::commonUtilities::setupEnvironment();
  flux_load();
  
  int count = default_count;
  std::string source_name(default_source);
  
  //TESTING MULTIPLE XML INPUT
  std::vector<std::string> fileList;
  fileList.push_back(facilities::commonUtilities::joinPath(facilities::commonUtilities::getXmlPath("GRBtemplate"), "GRBtemplate_library.xml"));
  FluxMgr fm(fileList);


  //FluxMgr fm;
  
  //Testing the addfactory function
  //    static PencilBeam* sean=TestSpec::instance();
  //    fm.addFactory("seantest", sean );
  //End Test
  
  if ( argn >1 ) source_name = argc[1];
  if( source_name =="help") { help(); return 0; }
  if( source_name =="list") { 
    listSources(fm.sourceList());
    listSpectra(); return 0; }
  if ( argn >2 ) count = ::atoi(argc[2]);
  
  cout << "------------------------------------------------------" << endl;
  cout << " Flux test program: type 'help' for help" << endl;
  cout << ( ( argn ==1)?  " No command line args, using default fluxes \""
	    :  " Selected source name \"");
  cout  << source_name <<"\"" << endl;
  
  
  std::list<std::string> source_list(fm.sourceList());
  
  if(( argn !=1) && std::find(source_list.begin(), source_list.end(), source_name)==source_list.end() ) {
    std::list<std::string> spectra(SpectrumFactoryTable::instance()->spectrumList());
    
    if( std::find(spectra.begin(), spectra.end(), source_name)==spectra.end() ) {
      std::cout << "Source \"" << source_name << "\" not found in the list or sources!" << std::endl;
      listSources(source_list);
      std::cout << "or in spectra list, which is:\n";
      listSpectra();
      
      return -1;
    }
  }
  std::list<std::string> allTheSources = fm.sourceList();
  std::list<std::string>::iterator abc;
  if(argn != 1){
    //    fm.test(std::cout, source_name, count);
    galacticTest(&fm,source_name);
    return 0;
  }
  std::string testfilename("test_GRB.out");
  std::ostream* m_out = new std::ofstream(testfilename.c_str());
  std::ostream& out = *m_out;
  std::cout << "Writing test results to the file " << testfilename << std::endl;
  int n=0;
  for(abc= allTheSources.begin() ; abc != allTheSources.end() && n++<5; abc++)
    {
      std::cout << "Source ["<<n<<"]: " << *abc << std::endl;
      out << "Source:  " << *abc <<std::endl;
      fm.test(out, (*abc), count);
      out << std::endl << std::endl << std::endl;
    }
  
  return 0;    
}

