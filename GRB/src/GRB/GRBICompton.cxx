#include "GRBICompton.h"
#include "GRBConstants.h"


GRBICompton::GRBICompton()
  : RadiationProcess()
{}

GRBICompton::GRBICompton(SpectObj spectrumObj)
  : RadiationProcess(spectrumObj)
{}


void GRBICompton::load(const GRBShock *Shock, 
		       const double time,
		       const double angle,
		       const double distance_to_source)
{
  load(time - Shock->tobs(),
       angle,
       distance_to_source,
       Shock->getGammaf(),
       Shock->getB(),
       Shock->getParticleN(),
       Shock->getGammaMin(),
       Shock->getGammaMax(),
       Shock->getVolume(),
       Shock->getThickness()
       );
}

void GRBICompton::load(const double time,
		       const double angle,
		       const double distance_to_source,
		       const double GAMMAF,
		       const double B,
		       const double N0,
		       const double gamma_min,
		       const double gamma_max,
		       const double Vol,
		       const double dr)
  
{
  if (time <= 0.){m_spectrumObj *= 0.; return;}
  
  int type = 0; // 0 -> Power Law Approximation

  double ComovingTime; 

  const double EnergyTransformation = (1.e+6*cst::mec2)*
    ( GAMMAF + sqrt(GAMMAF*GAMMAF+1.)*cos(angle) );
    
  //////////////////////////////////////////////////
  // Gamma e+-
  //////////////////////////////////////////////////
  
  std::vector<double>   x  = m_spectrumObj.getEnergyVector(1./(EnergyTransformation));
  std::vector<double>   dx = m_spectrumObj.getBinVector(1./(EnergyTransformation));
  std::vector<double> fsyn(cst::enstep+1, 0.0);
 
  const double Psyn_0    = 4./3.*Umag(B)*cst::erg2MeV*cst::c*cst::st; //MeV/s      
  const double esyn_0    = (3.*cst::pi/8.)*(B/cst::BQ)*cst::mec2; //MeV
  const double Psyn_esyn = Psyn_0/esyn_0; //1/sec
  const double tau_0     = (cst::flagIC == 1.) ? cst::c*cst::st*N0/Vol : cst::flagIC; //1/s
  
  if (type==0)
    {
      std::vector<double>::iterator x1 = x.begin();
      std::vector<double>::iterator its = fsyn.begin();
      while(x1!=x.end()) 
	{
	  ComovingTime = comovingTime(time,GAMMAF,(*x1),distance_to_source);
	  //comoving time <0 means that dispersive effect delayed 
	  //to much the photon at this energy bin
	  if(ComovingTime<0){its++; x1++; continue;}

	  //em and ec are specific: the compton scattering is with the gamma_min electron only...
	  double gc = cst::mec2/(Psyn_0*ComovingTime);
	  gc = (gc <= 1. ? 1.+1.e-6 : gc);
	  double em = 2.*pow(gamma_min,2.) *esyn_0*(pow(gamma_min,2.)-1.)/cst::mec2;
	  double ec = 2.*pow(gc,2.)*esyn_0 *(pow(gc,2.)-1.)/cst::mec2;
      
	  double gKN     = sqrt((cst::mec2/esyn_0)/(2.*gamma_min)+1.);	  
	  double gi      = sqrt(1. + (*x1/(2.*pow(gamma_min,2.)))*cst::mec2/esyn_0);
	  
	  double gi2     = pow(gi,2.);
	  double Psyn    = Psyn_0*(gi2-1.); //MeV/s
	  double tsyn    = ((gi)*cst::mec2)/Psyn;
	  double tsyn0    = (gamma_min*cst::mec2)/Psyn;
	  double tcross  = dr/cst::c;
	  double tau = tau_0;
	  tau *= (tsyn<tcross)? tsyn: tcross;
	  if(tau>1.) tau = 1.;
	  
	  double N_e;
	  if (gi > gKN)  
	    N_e = 0.0;
	  else
	    N_e = tau*N0;

	  
	  (*its) = processFlux(*x1,ec,em);
	  
	  if(ComovingTime <= tcross)
	    N_e *= tsyn0/tcross*(1.-exp(-ComovingTime/tsyn0));
	  else
	    N_e *= tsyn0/tcross*(1.-exp(-tcross/tsyn0))*exp((tcross-ComovingTime)/tsyn0);
	  (*its++) *= N_e; // adim
	  
	  (x1++);
	}
    }
  //////////////////////////////////////////////////
  // COMMON 
  //////////////////////////////////////////////////
  m_spectrumObj.clear();
  std::vector<double>::iterator x1  =    x.begin();
  std::vector<double>::iterator dx1 =   dx.begin();
  std::vector<double>::iterator its = fsyn.begin();
  while(dx1!=dx.end()) 
    {
      (*x1 ) *= EnergyTransformation; //eV
      (*dx1) *= 1.0e-6*EnergyTransformation; //MeV
      (*its) *= GAMMAF*(Psyn_esyn)/(*dx1); //1/s/MeV
      x1++;
      dx1++;
      its++;
    }
  (*x1) *= EnergyTransformation; //eV 
  (*its) = 0.0;
  m_spectrumObj.SetSpectrum(x,fsyn);
  return;
}



