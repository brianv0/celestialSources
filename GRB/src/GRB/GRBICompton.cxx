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
  //int type = 0; // 0 -> Power Law, from Sari and Esin, astro-ph/0005253
  //int type = 1; // 0 -> Power Law, the electron have gamma_min, Low IC
  //int type = 2; // 0 -> Power Law, The electron have gc, High IC
  int type = 3; // 0 -> Integration, Slow computation but all the 
  //                                 electrons are considered
  
  double ComovingTime; 
  
  const double EnergyTransformation = (1.e+6*cst::mec2)*
    ( GAMMAF + sqrt(GAMMAF*GAMMAF+1.)*cos(angle) );  // mc2->eV
  
  double   dgamma     = pow((gamma_max/gamma_min),1./cst::enstep);
  //////////////////////////////////////////////////
  // Gamma e+-
  //////////////////////////////////////////////////
  std::vector<double>   
    x    = m_spectrumObj.getEnergyVector(1./(EnergyTransformation)); //eV -> mc2
  std::vector<double>   
    dx   = m_spectrumObj.getBinVector(1.e+6/(EnergyTransformation));    //MeV -> mc2
  std::vector<double> 
    fsyn = m_spectrumObj.getSpectrumVector(); // is in ph/s/MeV
  std::vector<double> fIC(cst::enstep+1, 0.0);
  
  std::vector<double>::iterator x1  =    x.begin();
  std::vector<double>::iterator dx1 =   dx.begin();
  std::vector<double>::iterator isy = fsyn.begin();
  std::vector<double>::iterator its =  fIC.begin();
  const double Urad      = cst::alphae/cst::alphab*Umag(B);
  const double Psyn_0    = 4./3.*Umag(B)*cst::erg2MeV*cst::c*cst::st; //MeV/s      
  const double Pic_0     = 4./3.*  Urad *cst::erg2MeV*cst::c*cst::st; //MeV/s      
  const double esyn_0    = (3.*cst::pi/8.)*(B/cst::BQ)*cst::mec2; //MeV
  const double Psyn_esyn = Psyn_0/esyn_0; //1/sec
  const double Pic_eic = Pic_0/esyn_0; //1/sec

  const double tau_0     = (cst::flagIC == 1.) ? cst::c*cst::st*N0/Vol : cst::flagIC; //1/s
  const double tcross  = dr/cst::c;
  
  while(dx1!= dx.end()) 
    {
      (*isy) *= (*dx1)/(GAMMAF*(Psyn_esyn)*Vol);  // Now fsyn is in ph/cm^3
      isy++;
      dx1++;
    }


  if (type==0) //Power law from Sari and Esin, astro-ph/0005253
    {
      x1 = x.begin();
      its = fIC.begin();
      
      while(x1!=x.end()) 
	{
	  ComovingTime = comovingTime(time,GAMMAF,EnergyTransformation*(*x1),distance_to_source);
	  //comoving time <0 means that dispersive effect delayed 
	  //to much the photon at this energy bin
	  if(ComovingTime<0){its++; x1++; continue;}

	  //em and ec are specific: 
	  // the synchrotron cooling LF is:
	  double gc = cst::mec2/(Psyn_0*ComovingTime); 
	  gc = (gc <= 1. ? 1.+1.e-6 : gc);
	  // The IC energy that conrespond to gamma_min is (em_syn*gamma_min^2)
	  double em = 4./3.*pow(gamma_min,2.) *esyn_0*(pow(gamma_min,2.)-1.)/cst::mec2;
	  // The IC energy that conrespond to gc is (ec_syn*gc^2)
	  double ec = 4./3.*pow(gc,2.)        *esyn_0*(pow(gc,2.)-1.)/cst::mec2;
	  // The KN limit is when the electron give all its energy to the photon:
	  double ekn= gc; 
	  
	  // the gamma of the original syn electron is:
	  double gi      = gamma_min; //sqrt(1. + (3.*(*x1)*cst::mec2)/(4.*pow(gamma_min,2.)*esyn_0));
	  double gi2     = pow(gi,2.);
	  // The IC cooling time of the original electron is:
	  double tic    = (gi*cst::mec2)/(Pic_0*(gi2-1.));  
	  
	  double tau = tau_0;
	  tau *= (tic<tcross)? tic: tcross;
	  if(tau>1.) tau = 1.;
	  
	  double N_e = N0/Vol*(tau)*electronNumber(gi, gamma_min, gamma_max,
					    dr, ComovingTime, tic, N0);
	  
	  if(*x1<ekn)
	    (*its) = processFlux(*x1,ec,em);
	  else 
	    (*its) = 0.0;
	  
	  (*its++) *= N_e*gi*(dgamma-1.); // adim
	  
	  (x1++);
	}
    }
  else if (type==1) // The electron have gamma_min, Low IC
    {
      x1 = x.begin();
      its = fIC.begin();
      while(x1!=x.end()) 
	{
	  ComovingTime = comovingTime(time,GAMMAF,EnergyTransformation*(*x1),distance_to_source);
	  //comoving time <0 means that dispersive effect delayed 
	  //to much the photon at this energy bin
	  if(ComovingTime<0){its++; x1++; continue;}
	  
	  //em and ec are specific: the compton scattering is with the gamma_min electron only...
	  double gc = cst::mec2/(Psyn_0*ComovingTime);
	  gc = (gc <= 1. ? 1.+1.e-6 : gc);
	  
	  double em = 4./3.*pow(gamma_min,2.) *esyn_0 *(pow(gamma_min,2.)-1.)/cst::mec2;
	  double ec = 4./3.*pow(gamma_min,2.) *esyn_0 *(pow(gc,2.)-1.)/cst::mec2;
	  double ekn= gamma_min;
	  // gamma of the seed electron
	  double gi      = gamma_min; //sqrt(1. + (3.*(*x1)*cst::mec2)/(4.*pow(gamma_min,2.)*esyn_0));
	  double gi2     = pow(gi,2.);
	  
	  double tic    = (gi*cst::mec2)/(Pic_0*(gi2-1.));  
	  double tau = tau_0;
	  tau *= (tic<tcross)? tic : tcross;
	  if(tau>1.) tau = 1.;
	  
	  
	  double N_e = N0/Vol*(tau)*electronNumber(gi, gamma_min, gamma_max,
					    dr, ComovingTime, tic, N0);
	  if(*x1<ekn)
	    (*its) = processFlux(*x1,ec,em);
	  else 
	    (*its) = 0.0;
	  
	  (*its++) *= N_e*gi*(dgamma-1.); // adim

	  
	  (x1++);
	}
    }
  else if (type==2) // The electron have gc, High IC
    {
      x1 = x.begin();
      fIC.begin();
      while(x1!=x.end()) 
	{
	  ComovingTime = comovingTime(time,GAMMAF,EnergyTransformation*(*x1),distance_to_source);
	  //comoving time <0 means that dispersive effect delayed 
	  //to much the photon at this energy bin
	  if(ComovingTime<0){its++; x1++; continue;}
	  
	  //em and ec are specific: the compton scattering is with the gamma_min electron only...
	  double gc = cst::mec2/(Psyn_0*ComovingTime);
	  gc = (gc <= 1. ? 1.+1.e-6 : gc);
	  
	  
	  double em = 4./3.*pow(gc,2.) *esyn_0 *(pow(gamma_min,2.)-1.)/cst::mec2;
	  double ec = 4./3.*pow(gc,2.) *esyn_0 *(pow(gc,2.)-1.)/cst::mec2;
	  double ekn= gc;
	  double gi      = gc;//sqrt(1. + (3.*(*x1)*cst::mec2)/(4.*pow(gc,2.)*esyn_0));
	  double gi2     = pow(gi,2.);
	  double tic    = (gi*cst::mec2)/(Pic_0*(gi2-1.));  
	  double tau = tau_0;
	  tau *= (tic<tcross)? tic: tcross;
	  if(tau>1.) tau = 1.;
	  

	  double N_e = N0/Vol*(tau)*electronNumber(gc, gamma_min, gamma_max,
					    dr, ComovingTime, tic, N0);
	  
	  if(*x1<ekn)
	    (*its) = processFlux(*x1,ec,em);
	  else 
	    (*its) = 0.0;
	  
	  (*its++) *= N_e*gi*(dgamma-1.); // adim
	  
	  (x1++);
	}
    }
  else if (type==3) // Complete integration
    {
      x1  =    x.begin();
      //dx1 =   dx.begin();
      its =  fIC.begin();
      while(x1!=x.end())  // Output photon energy
	{ 
	  ComovingTime = comovingTime(time,GAMMAF,EnergyTransformation*(*x1),
				      distance_to_source);
	  if(ComovingTime<=0){its++; x1++; continue;}
	  for (int j = 0; j < cst::enstep;) // gamma electron
	    {
	      double gi      = gamma_min*pow(dgamma,j); 
	      double gi2     = pow(gi,2.);
	      //	      double temp    = (eic_0*(gi2-1.))/cst::mec2;
	      // temp = (temp>0.0) ? 1./temp : 0.0; 
	      double tic    = (gi*cst::mec2)/(Pic_0*(gi2-1.));  
	      double tau = tau_0; // adim
	      tau *= (tic<tcross)? tic: tcross;
	      if(tau>1.) tau = 1.;
	      double N_e = (tau)*electronNumber(gi, gamma_min, gamma_max,
						dr, ComovingTime, tic, N0);
	      int e0=0;
	      isy=fsyn.begin();
	      while(isy<fsyn.end()) // Initial energy
		{
		  (*its) += N_e*(*isy) 
		    *gamma_min*(pow(dgamma,j+2)-pow(dgamma,j))
		    *InverseComptonFunction(gi,x[e0],(*x1));
		  //cout<<" gamma  "
		  e0++;
		  isy+=2; // To speed up the integration. 
		}
	      j+=2;
	    }
	  its++;
	  x1++;
	}
    }
  //////////////////////////////////////////////////
  // COMMON 
  //////////////////////////////////////////////////
  m_spectrumObj.clear();
  x1  =    x.begin();
  dx1 =   dx.begin();
  its = fIC.begin();
  while(dx1!=dx.end()) 
    {
      (*x1 ) *= EnergyTransformation; //eV
      (*dx1) *= 1.0e-6*EnergyTransformation; //MeV
      (*its) *= GAMMAF*(Pic_eic)/(*dx1); //1/s/MeV
      
      x1++;
      dx1++;
      its++;
    }
  (*x1) *= EnergyTransformation; //eV 
  (*its) = 0.0;
  m_spectrumObj.SetSpectrum(x,fIC);
  return;
}

double GRBICompton::InverseComptonFunction(double g0,double e0,double energy)
{
  // note that e0 and energy are in unit of mec2
  double x0 = g0*e0;
  double q  = energy/(4.*g0*g0*e0)/(1.-(energy/g0));
  // ALSO E_MIN AND E_MAX
  double e_min= e0/(4.*g0*g0)/(1.+(e0/g0));
  double e_max= 4.*g0*g0*e0/(1.+(4.*g0*e0));
  
  if(( energy >= e_max) || (energy <= e_min)) return 0.;
  
  double f0 = (pow(4.*x0*q,2.)/(8.*x0*q+2.)+2.*q+1.)*(1.-q)+2.*q*log(q);
  //cout<<f0*energy/(4.*g0*g0*e0)<<endl;
  return f0*energy/(4.*g0*g0*e0);
}

