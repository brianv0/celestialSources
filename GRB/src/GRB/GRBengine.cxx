#include "GRBShell.h"
#include "GRBShock.h"
#include "GRBConstants.h"
#include "GRBengine.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"

GRBengine::GRBengine(GRBConstants *myParam)
{
  int type=3;
  switch(type)
    {
    case 1: // jet, no shell evolution
      {
	cout<<" CASE 1: jet like shells, no shell evolution, no collision, physical parameters "<<endl;  
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
	    double gamma_min = myParam->GammaMin();
	    //	    double gamma_max = myParam->GammaMax();
	    double gamma = gamma_min;
	    double mass = (1.-1.0e-3)*ei/(gamma*cst::c2); // Shell mass
	    
	    GRBShell ShockedMaterial(gamma,mass,myParam->Thickness(),myParam->JetRadius(),0.0,"jet");	  
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
	  }
      }
      break;
    case 2: // 2 shells that shock, jet, no shell evolution
      {
	cout<<"CASE 2: jet like shells, no shell evolution, collision , physical parameters "<<endl;
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
	    double gamma_min = myParam->GammaMin();
	    double gamma_max = myParam->GammaMax();
	    
	    double gamma_front = gamma_min+(gamma_max-gamma_min)*RandFlat::shoot(1.0);
	    double gamma_back  = gamma_front + (gamma_max-gamma_front)*RandFlat::shoot(1.0);
	    double mass_front = .5*ei/(gamma_front*cst::c2); // Shell mass
	    double mass_back =  .5*ei/(gamma_back*cst::c2); // Shell mass
	    GRBShell fShell(gamma_front,mass_front,myParam->Thickness(),myParam->JetRadius(),0.0,"jet");	  
	    GRBShell bShell(gamma_back ,mass_back ,myParam->Thickness(),myParam->JetRadius(),0.0,"jet");	  
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
      break;
    case 3: // jet, no shell evolution, obseravbles fixed
      {
	cout<<" CASE 3: jet like shells, no shell evolution, no collision, observable fixed "<<endl;  
	double rise_time   = myParam->RiseTime(); //sec
	double decay_time  = myParam->DecayTime(); //sec
	double peak_energy = myParam->PeakEnergy(); //MeV
	int NumberOfShocks       = myParam->Nshock();
	double ei = myParam->Etot()/NumberOfShocks; // internal energy of the shocked material
	//////////////////////////////////////////////////
	const double C1=3.3e-11;
	const double C2=cst::pi*1.e+10;
	const double C3=7.21787e-15;
	
	double gamma       = pow(pow(peak_energy,2.)*decay_time/(C3*C3*C2),1./5.);
	double thickness    = rise_time*gamma/C1;
	double jet_radius  = sqrt(ei*gamma*decay_time/(thickness*C2));
	myParam->setGammaMin(gamma);
	myParam->setGammaMax(gamma);
	myParam->setThickness(thickness);
	myParam->setJetRadius(jet_radius);
	
	m_duration=0.0;
	double BurstDuration     = myParam->Duration();
	if (BurstDuration<=0) BurstDuration=getDurationFromBATSE();
	double timeBetweenShocks = 2.*BurstDuration/NumberOfShocks*RandFlat::shoot(1.0);
	double ShockTime         = 0.;
	
	int ns=0;
	
	while (ns < NumberOfShocks)
	  {
	    double ComputedShockDuration = 0.0;
	    double mass = (1.-1.0e-3)*ei/(gamma*cst::c2); // Shell mass
	    //	    double internalEnergy  = ei/
	    GRBShell ShockedMaterial(gamma,mass,thickness,jet_radius,0.0,"jet");	  
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
	    ShockTime += timeBetweenShocks; 
	    timeBetweenShocks = 2.*BurstDuration/NumberOfShocks*RandFlat::shoot(1.0);
	    cout<<" **************************************************"<<endl;
	    cout<<" Estimated Rise Time      = "<<C1*thickness/gamma<<endl;
	    cout<<" Estimated Decay Time     = "<<C2*pow(jet_radius,2.)*thickness/(gamma*ei)
		<<endl;
	    cout<<" Estimated Break   Energy   = "<<C3*pow(gamma,3.)*sqrt(ei/thickness)/(jet_radius)<<endl;
	  }
      }
      break;
    }
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
