#include "GRBShell.h"
#include "GRBShock.h"
#include "GRBConstants.h"
#include "GRBengine.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"

GRBengine::GRBengine(GRBConstants *myParam)
{

  //////////////////////////////////////////////////
  const double C1=3.3e-11*1.01688;
  const double C2=cst::pi*1.e+10*46.0;
  const double C3=2.90851206811141e-14/0.3
    *pow((cst::p-2.)/(cst::p-1.),2.)/(2./(1+cos(myParam->JetAngle()))); 
  
  //////////////////////////////////////////////////

  int engine_type = myParam->EngineType();

  std::string shell_type = (myParam->ShellType()==1) ? "jet" : "iso";
  double gamma_min = myParam->GammaMin();
  double gamma_max = myParam->GammaMax();
  
  if(engine_type==1)
    {
      std::cout<<
	" CASE 1: no shell evolution, no collision, physical parameters "
	       <<std::endl;  
      m_duration=0.0;
      double BurstDuration     = myParam->Duration();
      if (BurstDuration<=0) BurstDuration=getDurationFromBATSE();
      int NumberOfShocks       = myParam->Nshock();
      double timeBetweenShocks = 
	2.*BurstDuration/NumberOfShocks*RandFlat::shoot(1.0);
      double ShockTime         = 0.;
      // internal energy of the shocked material
      double ei = myParam->Etot()/NumberOfShocks; 
      int ns=0;
      
      while (ns < NumberOfShocks)
	{
	  double  ShockDuration = 0.1;
	  if(NumberOfShocks==1) ShockDuration = BurstDuration;
	  double ComputedShockDuration = 0.0;
	  double gamma = gamma_min+(gamma_max-gamma_min)*RandFlat::shoot(1.0);
	  double mass = (1.-1.0e-3)*ei/(gamma*cst::c2); // Shell mass
	  
	  GRBShell ShockedMaterial(gamma,mass,
				   myParam->Thickness(),
				   myParam->JetRadius(),
				   myParam->ShellRadius(),
				   shell_type);	  
	  ShockedMaterial.setEint(1.0e-3*ei);
	  GRBShock iShock(ShockedMaterial);
	  
	  ComputedShockDuration = iShock.duration();
	  iShock.setTobs(ShockTime);
	  theShocks.push_back(iShock);		
	  if (ShockTime+ComputedShockDuration>m_duration)
	    {
	      m_duration = (ShockTime+ComputedShockDuration < 1.5*BurstDuration) ? 
		ShockTime+ComputedShockDuration : 1.5*BurstDuration;
	    }
	  ns++;
	  ShockTime += timeBetweenShocks; 
	  timeBetweenShocks = 2.*BurstDuration/NumberOfShocks*RandFlat::shoot(1.0);
	  std::cout<<" **************************************************"<<std::endl;
	  std::cout<<" Estimated Rise Time      = "<<
	    C1*myParam->Thickness()/gamma<<std::endl;
	  std::cout<<" Estimated Decay Time     = "<<
	    C2*pow(myParam->JetRadius(),2.)*myParam->Thickness()/(gamma*gamma*ei)<<std::endl;
	  std::cout<<" Estimated Break Energy   = "<<
	    C3*pow(gamma,3.)*sqrt(ei/myParam->Thickness())/(myParam->JetRadius())<<std::endl;
	}
    }
  else if(engine_type == 2) 
    {
      std::cout<<
	"CASE 2: no shell evolution, collision , physical parameters "
	       <<std::endl;
      m_duration=0.0;
      double BurstDuration     = myParam->Duration();
      if (BurstDuration<=0) BurstDuration=getDurationFromBATSE();
      int NumberOfShocks       = myParam->Nshock();
      double timeBetweenShocks = 2.*BurstDuration/NumberOfShocks*RandFlat::shoot(1.0);
      double ShockTime         = 0.;
      double ei = myParam->Etot()/NumberOfShocks; // internal energy of the shocked material
      int ns=0;
      
      while (ns < NumberOfShocks)
	{
	  double  ShockDuration = 0.1;
	  if(NumberOfShocks==1) ShockDuration = BurstDuration;
	  double ComputedShockDuration = 0.0;
	  
	  double gamma_front = gamma_min+(gamma_max-gamma_min)*RandFlat::shoot(1.0);
	  double gamma_back  = gamma_front + (gamma_max-gamma_front)*RandFlat::shoot(1.0);
	  double mass_front  = .5*ei/(gamma_front*cst::c2); // Front Shell mass
	  double mass_back   = .5*ei/(gamma_back *cst::c2); // Back  Shell mass
	  GRBShell fShell(gamma_front,mass_front,
			  myParam->Thickness(),
			  myParam->JetRadius(),
			  myParam->ShellRadius(),
			  shell_type);	  
	  GRBShell bShell(gamma_back ,mass_back,
			  myParam->Thickness(),
			  myParam->JetRadius(),
			  myParam->ShellRadius(),
			  shell_type);	  
	  
	  GRBShell ShockedMaterial = fShell + bShell;
	  GRBShock iShock(ShockedMaterial);
	  ComputedShockDuration = iShock.duration();
	  iShock.setTobs(ShockTime);
	  theShocks.push_back(iShock);		
	  
	  if (ShockTime+ComputedShockDuration>m_duration)
	    {
	      m_duration = (ShockTime+ComputedShockDuration < 1.5*BurstDuration) ? 
		ShockTime+ComputedShockDuration : 1.5*BurstDuration;
	      
	    }
	  ns++;
	  ShockTime += timeBetweenShocks; 
	  timeBetweenShocks = 2.*BurstDuration/NumberOfShocks*RandFlat::shoot(1.0);
	  
	  
	  
	}
    }
  else if(engine_type == 0) 
    {
      std::cout<<
	" CASE 0: no shell evolution, no collision, observable fixed "
	       <<std::endl;  
      int NumberOfShocks       = myParam->Nshock();
      m_duration=0.0;
      double BurstDuration     = myParam->Duration();
      if (BurstDuration<=0) BurstDuration=getDurationFromBATSE();
      
      double ShockTime         = 0.;
      
      int ns=0;
      
      while (ns < NumberOfShocks)
	{
	  double rise_time   = (2.*RandFlat::shoot(1.0))* myParam->RiseTime(); //sec
	  double decay_time  = (2.*RandFlat::shoot(1.0))* myParam->DecayTime();//myParam->DecayTime(); //sec
	  double peak_energy = myParam->PeakEnergy(); //MeV
	  
	  double ei = (2.*RandFlat::shoot(1.0))*myParam->Etot()/NumberOfShocks; // internal energy of the shocked material
	  
	  double gamma             = pow(pow(peak_energy,2.)*decay_time/(C3*C3*C2),1./4.);
	  double thickness         = rise_time*gamma/C1;
	  double jet_radius        = gamma*sqrt(ei*decay_time/(thickness*C2));
	  double shell_radius      = 0.5*jet_radius;
	  
	  myParam->setGammaMin(gamma);
	  myParam->setGammaMax(gamma);
	  myParam->setThickness(thickness);
	  myParam->setJetRadius(jet_radius);
	  myParam->setShellRadius(shell_radius);
	  
	  
	  double ComputedShockDuration = 0.0;
	  double mass = (1.-1.0e-3)*ei/(gamma*cst::c2); // Shell mass
	  //	    double internalEnergy  = ei/
	  GRBShell ShockedMaterial(gamma,mass,
				   thickness,
				   jet_radius,
				   shell_radius,
				   shell_type);	  
	  ShockedMaterial.setEint(1.0e-3*ei);
	  GRBShock iShock(ShockedMaterial);
	  
	  ComputedShockDuration = iShock.duration();
	  iShock.setTobs(ShockTime);
	  theShocks.push_back(iShock);		
	  if (ShockTime+ComputedShockDuration>m_duration)
	    {
	      m_duration = ShockTime+ComputedShockDuration;
	    }
	  ns++;
	  double timeBetweenShocks = 2.*BurstDuration/NumberOfShocks*RandFlat::shoot(1.0);
	  ShockTime += timeBetweenShocks; 
	  
	  std::cout<<" **************************************************"<<std::endl;
	  std::cout<<" Estimated Rise Time      = "<<C1*thickness/gamma<<std::endl;
	  std::cout<<" Estimated Decay Time     = "<<C2*pow(jet_radius,2.)*thickness/(gamma*gamma*ei)<<std::endl;
	  std::cout<<" Estimated Break Energy   = "<<C3*pow(gamma,3.)*sqrt(ei/thickness)/(jet_radius)<<std::endl;
	}
    }
  else
    {
      std::cout<<" No valid option in GRBParam engine type."<<std::endl;
      std::cout<<" Possible choices are:"<<std::endl;
      std::cout<<" 0 Physical Parameter, for the shocked matherial "
	       <<std::endl;
      std::cout<<" 1 Observed Parameter fixed "<<std::endl;
      std::cout<<
	" 2 Physical Parameter, for the shells, Shock from their collision "
	       <<std::endl;
      exit(0);
    }
  std::cout<<""<<std::endl;
  if(cst::flagQG) std::cout<<" Quantum Gravity is ON !!"<<std::endl;
  if(cst::flagIC>0) std::cout<<" Inverse Compton is ON !!"<<std::endl;
  std::cout<<""<<std::endl;

  //Glast Direction
  //  m_direction=std::make_pair(((RandFlat::shoot(1.0))*1.4)
  //  			  -0.4,(RandFlat::shoot(1.0))*2*M_PI);
  // Galactic (l,b)
  m_direction = std::make_pair(((RandFlat::shoot(1.0))*360.0)-180.,
			       (RandFlat::shoot(1.0)*180.0)-90.);
  m_distance  = myParam->Distance();
}

double GRBengine::getDurationFromBATSE(char* burst_type)
{
  double dur;
  //////////////////////////////////////////////////
  // It determines if, in case of random selection of the parameters,
  // the burst is long or short...
  if (burst_type!="Short" && burst_type!="Long")
    {
      if (RandFlat::shoot(1.0)<=0.3)
	{burst_type="Short";}
      else
	{burst_type="Long";}
    }  
  if(burst_type=="Short")
    {
      double temp=RandGauss::shoot(-3.65755e-01,5.14102e-01);//GRBConstants::SelectGaussRandom(0.0,2.5);
      dur=pow(10.,temp);
    }
  else
    {
      double temp=RandGauss::shoot(1.46819e+00,4.91505e-01);//GRBConstants::SelectGaussRandom(1.0,3.0);
      dur=pow(10.,temp);
    }
  return dur;
}


