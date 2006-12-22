/*!\file GRBTest.cxx
 * \brief test program for GRB simulation studies.
 * 
 * This executable uses ROOT to display several histograms at run time.
 * To use it, type on the prompt ./GRBTest.exe
 * After some output, an empty canvas named \b c1 will appear, 
 * and the program hangs at the line:
 * "Type a time (in sec) or 0 for the complete evolution"
 * \arg A time \c t entered will provoke the computation of the flux at \c t
 * seconds after initiation of the burst in the GLAST referential frame.
 * Then, the canvas \c1 will show on top the number of photons reaching the 
 * detector as a function of their energy (in log scale eV), and on bottom the
 * sampled photons from the preceding curve (note that the lowest energy for 
 * sampling is set to 10 Mev.)
 * \arg If 0 is typed, the complete simulation of GRB photon arrival to the 
 * detector is simulated. The canvas displays the continuous evolution of 
 * the spectrum. At the end, 4 new canvas pops up, showing several summary 
 * histograms.
 */
#include <iterator>
#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>
//Include files for spectrum...
#include "src/FluxMgr.h"
#include "FluxSvc/EventSource.h"
#include "src/SpectrumFactoryTable.h"
//include files for ROOT...
#include "TTree.h"
#include "TBranch.h"
#include "TObjArray.h"
#include "TFile.h"

//GRB include files...
#include "../../src/GRB/GRBShell.h"
#include "../../src/GRB/GRBShock.h"
#include "../../src/GRB/GRBConstants.h"
#include "../../src/GRB/GRBSim.h"
#include "../../src/GRB/GRBSpectrum.h"

using namespace std;

//! The size of the energy array is taken from GRBConstant.h
#define ENERGYSIZE cst::enstep
//! The size of the time array is taken from GRBConstant.h
#define TIMESIZE cst::nstep
/*! the flux is a bidimentional array (matrix?) funtion of energy and time.
 * Its Unities are 
 *
 */
const double TIME=10.0;
const double EVENTS=100000;
double m_flux[ENERGYSIZE][TIMESIZE];
double timex[TIMESIZE];
double energyx[ENERGYSIZE];

static const char * default_arg="GRBSpectrum";
int test1(int argc, char** argv);

void help() {
   std::cout << 
      "   Simple test program for Transient Sources.\n"
      "   Command line args are \n"
      "      '-events <number of events to create>'\n"
      "      '-time <time in seconds>'    for the maximum time\n"
      "      '-list' lists the available spectra\n"
      "      '-help' for this help"
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

#define DLL_DECL_SPECTRUM(x)   extern const ISpectrumFactory& x##Factory; x##Factory.addRef();


void flux_load() {
  // these are the spectra that we want to make available
  DLL_DECL_SPECTRUM( CHIMESpectrum);
  DLL_DECL_SPECTRUM( AlbedoPSpectrum);
  DLL_DECL_SPECTRUM( HeSpectrum);
  DLL_DECL_SPECTRUM( GalElSpectrum);   
  //   DLL_DECL_SPECTRUM( CrElectron);
  //DLL_DECL_SPECTRUM( CrProton);
  DLL_DECL_SPECTRUM( FILESpectrum);
}


int main(int argc, char** argv)
{
  test1(argc,argv);
  return 0;
}

double CalculateFluence(double ee,double e1=cst::enmin*1e-6,double e2=cst::enmax*1e-6)
{
  double fluence;
  if (ee>=e1 && ee<=e2) 
    {
      fluence=ee;
    }
  else 
    {
      fluence=0.0;
    }
  return fluence;
}

int test1(int argc, char** argv)
{
  int nume,i;
  double dt,DeadTime;
  int num_sources=0;
  double time_max=TIME;  //time to use for flux and rate functions
  double events_max=EVENTS;
  double time,energy,Rate;
  Vector dir;
  double cos_theta,phi;
  int current_arg = 1;
  double fluence1;
  std::string arg_name(default_arg);
  vector<std::string> sources;
  
  //flux_load();
  
  FluxMgr fm(sources); 
  
  cout << "------------------------------------------------------" <<endl;
  cout << " Flux test program: type 'GRBTest.exe -help' for help" <<endl;
  cout << ( ( argc == 1)?  " No command line args, using defaults"
      :  "") <<endl;

  while(current_arg < argc)
    {
      arg_name = argv[current_arg];
      if("-help" == arg_name || "help" == arg_name) 
	{ 
	  help();
	  return 0;
	}
      
      else if("-list" == arg_name) 
	{ 
	  listSources(fm.sourceList());
	  listSpectra(); 
         return 0; 
      }

      else if("-time" == arg_name) {
	time_max = atof(argv[++current_arg]);
	cout<<" MAX TIME = "<<time_max<<endl;
      }
      else if("-events" == arg_name) {
	events_max = atof(argv[++current_arg]);
	cout<<" MAX NUM OF EVENTS = "<<events_max<<endl;
	if (events_max<=1) return 0;
      }
      else if('-' == arg_name[0]) {std::cerr << "Unrecognized option "<< arg_name << ", -help for help" << std::endl;}
      else
	{
	  sources.push_back(arg_name);
	  num_sources++;
	}
      current_arg++;
   }
  if(0 == sources.size())
    {
      sources.push_back(default_arg);
      num_sources++;
    }
  cout<<"Num Sources = "<<num_sources<<endl;

  // Create the file, the tree and the branches...
  TObjArray Forest(0);
  TTree* events;        
  
  const char* name;
  // Here you can set the DeadTime of the instrument:
  DeadTime=1.0e-5; //Sec
  //
  for(i = 0; i < num_sources; i++)
    {
      nume=0;
      fluence1=0.0;
      EventSource *e = fm.source(sources[i]);
      cout<<" Source Name = "<<sources[i]<<endl;
      name=sources[i].c_str();
      
      events= new TTree(name,name);
      events->Branch("energy",&energy,"energy/D");
      events->Branch("time",&time,"time/D");
      events->Branch("Rate",&Rate,"Rate/D");
      events->Branch("cos_theta",&cos_theta,"cos_theta/D");
      events->Branch("phi",&phi,"phi/D");
      Forest.Add(events);
      

      pair<double,double> loc=fm.location();
      cout << loc.first << "   " << loc.second <<endl;
      time=1.0e-4;

      while (time<=time_max && nume<=events_max)
	{
	  FluxSource *f = e->event(time);
	  dir = f->launchDir();
	  
	  cos_theta = dir.z();
	  phi = atan2(dir.y(),dir.x());
	  if(phi < 0) phi = 2*M_PI + phi;

	  phi=180*phi/M_PI; //degrees

	  Rate= e->rate(time);
	  dt=(e->interval(time))/(e->totalArea());
	  if (dt<DeadTime) dt=DeadTime;
	  energy = f->energy();
	  fluence1+=CalculateFluence(energy);
	  events->Fill();
	  if (nume%1==0){
	    cout<<
	      "-------- Event Number: "<<nume<<"\n"<<
	      " Time [s] = "<<time<<"\n"<<
	      " Rate [ph/(s)]= "<<Rate<<"\n"<<
	      // " Flux [ph/(m^2 s)]= "<<e->flux(time)<<"\n"<<
	      " 1/Rate [s]= "<<1/Rate<<"\n"<<
	      " Interval [s]= "<<dt<<"\n"<<
	      " Area [m^2]= "<<e->totalArea()<<"\n"<<
	      // " --------------------------------\n"<<
	      " Energy of Photon Extracted [ MeV ]= "<<energy<<"\n"<<
	      " Direction: Cos(theta) = " << cos_theta <<", phi = "<<phi<<"\n"<<
	      endl;
	  }
	  time+=dt;
	  nume++;
	}
	cout<<"Time final="<<time<<endl;
	cout<<"Number of events processed ="<<nume<<endl;
	events->Print();
	delete e;
	cout<<"Fluence [erg/cm^2]="<<fluence1/(e->totalArea()*1.0e+4)*(1.0/cst::erg2MeV)<<endl;
	//	events->Scan("energy:time:Rate:cos_theta:phi");
    } // For all the sources...
  TFile f("Events.root","recreate");
  Forest.Write();
  f.Close();
}
