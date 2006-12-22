//#include <iostream.h>
//#include <fstream.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <math.h>
#include "src/FluxException.h" // defines FATAL_MACRO
#include "FluxSvc/mainpage.h"
#include "GRBConstants.h"

#include "facilities/Util.h"

using namespace std;

GRBConstants::GRBConstants()
{
  ReadParam();
}

void GRBConstants::ReadParam(){
  cout<<"read parameters"<<endl;
  char buf[100];
  double g0,g1;
  std::string paramFile = "$(GRBROOT)/src/test/GRBParam.txt";
  facilities::Util::expandEnvVar(&paramFile);
  ifstream f1(paramFile.c_str());
  if (! f1.is_open()) 
    {
      cout<<"Error Opening $(GRBROOT)/src/test/GRBParam.txt\n";
      //TODO LIST: still need to remove this exit, without gwtting a core dump!
      exit(1);
    }
  cout<<"Read the file: "<<paramFile.c_str()<<endl;
  
  f1.getline(buf,100);
  sscanf(buf,"%d",&nshell);
  setNshell(nshell);
  
  f1.getline(buf,100);
  sscanf(buf,"%lf",&redshift);
  setRedshift(redshift);
  
  f1.getline(buf,100);
  sscanf(buf,"%lf",&etot);
  setEtot(etot);
  
  f1.getline(buf,100);
  sscanf(buf,"%lf",&r0);
  setR0(r0);
  
  f1.getline(buf,100);
  sscanf(buf,"%lf",&t0);
  setT0(t0);
  
  f1.getline(buf,100);
  sscanf(buf,"%lf",&g0);
  setGamma0(g0);
  
  f1.getline(buf,100);
  sscanf(buf,"%lf",&g1);
  setDGamma(g1);
  cout<<g0<<"  "<<g1<<endl;
  
      
  f1.close();

}

