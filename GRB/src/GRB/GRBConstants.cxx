#include <iostream>
#include <stdio.h>
#include <string>
#include <ctime>

#include "GRBConstants.h"
#include "facilities/Util.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RanluxEngine.h"

using namespace std;

GRBConstants::GRBConstants()
{
  ReadParam();
}

// To initialize the random number generator with a different seed each run, we use the clock time.
// The seed can be also initialized to a given number.
void GRBConstants::InitializeRandom(long seed)
{
  HepRandom::setTheEngine(new RanluxEngine);
  if(seed == 0)
    {
      time_t ctime;
      if (time(&ctime)==time_t(-1))
	{
	  std::cout<<" Time() not defined for this processor "<<std::endl;
	  seed=1;
	}
      seed=(ctime);
    }
  HepRandom::setTheSeed(seed);
  std::cout<<"HepRandom initialized... SEED = "<<seed<<std::endl;
  HepRandom::showEngineStatus();
}

HepRandomEngine* GRBConstants::GetTheRandomEngine(long seed)
{
  InitializeRandom(seed);
  return HepRandom::getTheEngine();
}
// This method read the parameter file, and sets all the private data members.
int GRBConstants::ReadParam(){   
  char buf[100];
  std::string paramFile = "$(GRBROOT)/src/test/GRBParam.txt";
  facilities::Util::expandEnvVar(&paramFile);
  ifstream f1(paramFile.c_str());
  if (! f1.is_open()) 
    {
      std::cout<<"Error Opening $(GRBROOT)/src/test/GRBParam.txt\n";
      exit(1);
    }

  f1.getline(buf,100);
  sscanf(buf,"%d",&m_nshock);

  f1.getline(buf,100);
  sscanf(buf,"%lf",&m_etot);

  f1.getline(buf,100);
  sscanf(buf,"%lf",&m_redshift);
  
  f1.getline(buf,100);
  sscanf(buf,"%lf",&m_enph);

  f1.getline(buf,100);
  sscanf(buf,"Shell type = %d",&m_shellType);
  
  f1.getline(buf,100); //## If Shell == 1 : Jet Shells
  
  f1.getline(buf,100);
  sscanf(buf,"%lf",&m_r0);
  
  f1.getline(buf,100);
  sscanf(buf,"%lf",&m_angle);

  f1.getline(buf,100); //## If Shell == 0 : Isotropic Shells
  
  f1.getline(buf,100);
  sscanf(buf,"%lf",&m_d0);

  f1.getline(buf,100);
  sscanf(buf,"Engine type = %d",&m_engineType); 
  
  f1.getline(buf,100); // ## If  Engine type == 0 (Observables pparameters are):
  
  f1.getline(buf,100);  
  sscanf(buf,"%lf",&m_duration);
  
  f1.getline(buf,100);
  sscanf(buf,"%lf",&m_rt);
  
  f1.getline(buf,100);
  sscanf(buf,"%lf",&m_dt);
  
  f1.getline(buf,100);
  sscanf(buf,"%lf",&m_peak);
  
  f1.getline(buf,100); //## If  Engine type == 1 or 2 (Observables pparameters are):

  f1.getline(buf,100);
  sscanf(buf,"%lf",&m_t0);
    
  f1.getline(buf,100);
  sscanf(buf,"%lf",&m_g0);   
  
  f1.getline(buf,100);
  sscanf(buf,"%lf",&m_g1);
  
  //////////////////////////////////////////////////
  setNshock(m_nshock);
  setEtot(m_etot);
  setRedshift(m_redshift);
  setEnergyPh(m_enph);
  
  //## If Shell == 1 : Jet Shells
  setShellType(m_shellType);
  setJetRadius(m_r0);
  setJetAngle(m_angle);

  //## If Shell == 0 : Isotropic Shells
  setShellRadius(m_d0);

  setEngineType(m_engineType);
  // ## If  Engine type == 0 (Observables pparameters are):
  setDuration(m_duration);
  setRiseTime(m_rt);
  setDecayTime(m_dt);
  setPeakEnergy(m_peak);
  //  ## If  Engine type == 1 or 2 (Observables pparameters are)
  setThickness(m_t0);
  setGammaMin(m_g0);
  setGammaMax(m_g1);

  return 0;
}

double GRBConstants::Distance()
{
  double qo=(1.0+3.0*cst::wzel)/2.0;
  return ((cst::c/(cst::Hubble*1.0e+5)/pow(qo,2.0))*
	  (m_redshift*qo+(qo-1.0)*(-1.0+sqrt(2.0*qo*m_redshift+1.0)))*cst::mpc2cm);
}

void GRBConstants::Print()
{
  // std::cout<<"Read the file: "<<paramFile.c_str()<<std::endl;
  std::cout<<"********** Selected Parameters: **********"<<std::endl;
  std::cout<<" Number of Shocks               = "<<Nshock()<<std::endl;
  std::cout<<" Total Energy at the Source     = "<<Etot()<<std::endl;
  std::cout<<" Redshift of the Source         = "<<Redshift()<<std::endl;
  std::cout<<" Minimum energy extracted       = "<<EnergyPh()<<std::endl;
  if(m_shellType==1)
    {
      std::cout<<""<<std::endl;
      std::cout<<" Shell type selected : JET "<<std::endl;
      std::cout<<"    Angle  of the Jet     = "<<JetAngle()<<std::endl;
    }
  if(m_shellType==0)
    {
      std::cout<<""<<std::endl;
      std::cout<<" Shell type selected : ISO "<<std::endl;
    }
  if(m_engineType==0)
    {
      std::cout<<""<<std::endl;
      std::cout<<" Engine type selected : Observed Parameters "<<std::endl;
      std::cout<<"    Duration of the Burst = "<<Duration()<<std::endl;
      std::cout<<"    Rise  time of a peak  = "<<RiseTime()<<std::endl;
      std::cout<<"    Decay time of a peak  = "<<DecayTime()<<std::endl;
      std::cout<<"    Peak  energy (MeV)    = "<<PeakEnergy()<<std::endl;
    } 
  if((m_engineType==1) || (m_engineType==2))
    {
      std::cout<<" Engine type selected : Physical Parameters "<<std::endl;
      if(m_shellType==1)     
	{  
	  std::cout<<"    Radius of the Jet(cm) = "<<JetRadius()<<std::endl;
	}
      if(m_shellType==0)
	{
	  std::cout<<"    Shell radius (cm)     = "<<ShellRadius()<<std::endl;
	}     
      std::cout<<"    Shell Thickness  (cm)  = "<<Thickness()<<std::endl;
      std::cout<<"    Minimum Lorentz factor = "<<GammaMin()<<std::endl;
      std::cout<<"    Maximum Lorentz factor = "<<GammaMax()<<std::endl;
    }
  std::cout<<"*******************************************"<<std::endl;
}

void GRBConstants::Save(bool flag)
{
  if(flag)
    {
      m_paramFile=cst::paramFile;
      facilities::Util::expandEnvVar(&m_paramFile);
      std::cout<<"Save Parameters into the file: "<<m_paramFile.c_str()<<std::endl;
      ofstream f2(m_paramFile.c_str(),ios::app);
      if (! f2.is_open()) 
	{
	  std::cout<<"Error Opening "<<m_paramFile<<std::endl;
	  exit(1);
	}
      
      f2<<Nshock()<<std::endl;
      f2<<Redshift()<<std::endl;
      f2<<Etot()<<std::endl;
      f2<<JetRadius()<<std::endl;
      f2<<Thickness()<<std::endl;
      f2<<GammaMin()<<std::endl;
      f2<<GammaMax()<<std::endl;
      f2.close();
    }
}

