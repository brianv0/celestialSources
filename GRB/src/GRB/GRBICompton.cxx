#include "GRBICompton.h"
#include "GRBConstants.h"


GRBICompton::GRBICompton(SpectObj spectrumObj)
  : m_spectrumObj(spectrumObj)
{}
void GRBICompton::load(const GRBShock *Shock, 
			  const double time,
			  const double angle)
{
  load(time - Shock->tobs(),
       angle,
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
			  const double GAMMAF,
			  const double B,
			  const double N0,
			  const double gamma_min,
			  const double gamma_max,
			  const double Vol,
			  const double dr)
  
{
  //  int type = 1; // integration
  int type = 0; // Power Law Approximation
  if (time <= 0.){m_spectrumObj*=0.; return;}
  {
    const int n = cst::enstep; 
    const double ComovingTime = time*GAMMAF; 
    const double EnergyTransformation = (1.e+6*cst::mec2)*(GAMMAF+sqrt(GAMMAF*GAMMAF+1.)*cos(angle));
    
    x.clear();
    g.clear();
    dx.clear();
    
    //////////////////////////////////////////////////
    // Gamma e+-
    //////////////////////////////////////////////////
    
    x  = m_spectrumObj.getEnergyVector(1./(EnergyTransformation));
    dx = m_spectrumObj.getBinVector(1./(EnergyTransformation));
    
    const double Psyn_0 = 4./3.*Umag(B)*cst::erg2MeV*cst::c*cst::st; //MeV/s      
    const double esyn_0 = (3.*cst::pi/8.)*(B/cst::BQ)*cst::mec2; //MeV
    const double Psyn_esyn = Psyn_0/esyn_0; //1/sec
    const double N_0  = N0*(cst::p-1.)/(pow(gamma_min,1.0-cst::p)-pow(gamma_max,1.0-cst::p)); 
    const double tau_0 = (cst::flagIC == 1.) ? cst::c*cst::st*N0/Vol : cst::flagIC; //1/s
    
    if (type==0)
      {
	double em = 2.*pow(gamma_min,2.) *esyn_0*(pow(gamma_min,2.)-1.)/cst::mec2;
	double gc = cst::mec2/(Psyn_0*ComovingTime);
	gc = (gc <= 1. ? 1.+1.e-6 : gc);
	double ec = 2.*pow(gc,2.)*esyn_0 *(pow(gc,2.)-1.)/cst::mec2;
	//cout<<" Ec IC (MeV) = "<<ec*(cst::mec2)<<" Em IC (MeV) = "<<em*(cst::mec2)<<endl;	 
	fsyn.clear();
	fsyn.resize(n+1,0.0);
	
	std::vector<double>::iterator x1 = x.begin();
	std::vector<double>::iterator its = fsyn.begin();
	while(x1!=x.end()) 
	  {
	    double gKN     = sqrt((cst::mec2/esyn_0)/(2.*pow(gamma_min,2.))+1.);
	    
	    double gi      = sqrt(1. + (*x1/(2.*pow(gamma_min,2.)))*cst::mec2/esyn_0);
	    
	    
	    double gi2     = pow(gi,2.);
	    double Psyn    = Psyn_0*(gi2-1.); //MeV/s
	    double tsyn    = ((gi)*cst::mec2)/Psyn;
	    double tcross  = dr/(sqrt(1.-1./gi2)*cst::c);
	    double tau = (tsyn<tcross)?tsyn*tau_0:tcross*tau_0;
	    tau = (tau<1.) ? tau:1.;
	    
	    double N_e = (tau)*N_0;
	    
	    if (gi<gamma_min ) 
	      {
		N_e *= gi/gamma_min*pow(gamma_min,-cst::p); //adim;
	      } 
	    else if (gi > gamma_max) 
	      {
		N_e *= 0.0;
	      }
	    else 
	      {
		N_e *= pow(gi,-cst::p); //adim;
	      }	      
	    if (gi > gKN)  N_e *= 0.0;
	    if(ComovingTime<tcross)
	      {
		//	N_e *= ComovingTime/tcross*exp(-ComovingTime/tsyn);
		N_e *= tsyn/tcross*(1.-exp(-ComovingTime/tsyn));
	      }
	    else
	      {
		N_e *= tsyn/tcross*(1.-exp(-tcross/tsyn))*exp((tcross-ComovingTime)/tsyn);
		//	N_e *= exp(-ComovingTime/tsyn);
	      }
	    
	    if(ec <= em) //FAST COOLING
	      {
		if((*x1) <= ec)
		  {
		    (*its) = (pow(*x1/ec,0.33333));
		  }
		else if ((*x1) <= em)
		  {
		    (*its) = (pow(*x1/ec,-0.5));
		  }
		else 
		  {
		    (*its) = (pow(em/ec,-0.5))*(pow(*x1/em,-cst::p/2.));
		  }
	      }
	    else //SLOW COOLING
	      {
		//		  cout<<"*";
		if((*x1) <= em)
		  {
		    (*its) = (pow(*x1/em,0.33333));
		  }
		else if ((*x1) <= ec)
		  {
		    (*its) = (pow(*x1/em,-(cst::p-1.)/2.));
		  }
		else 
		  {
		    (*its) = (pow(ec/em,-(cst::p-1.)/2.)*(pow(*x1/ec,-cst::p/2.)));
		  }
	      }
	    (*its++) *= N_e; // adim
	    (x1++);
	  }
      }
    else
      {
	//////////////////////////////////////////////////
	// Integration
	////////////////////////////////////////////////// 	  
	
	fsyn.clear();
	fsyn.resize(n+1,0.0);
	std::vector<double>::iterator x1 = x.begin();
	std::vector<double>::iterator its = fsyn.begin();
	
	while(x1!=x.end()) 
	  {
	    double gi      = sqrt(1. + (*x1)*cst::mec2/esyn_0);
	    double gi2     = pow(gi,2.);
	    double Psyn    = Psyn_0*(gi2-1.); //MeV/s
	    double esyn    = esyn_0*(gi2-1.); //MeV
	    double tsyn    = ((gi)*cst::mec2)/Psyn;
	    double tcross  = dr/(sqrt(1.-1./gi2)*cst::c);
	    double tau = (tsyn<tcross)?tsyn*tau_0:tcross*tau_0;
	    tau = (tau<1.) ? tau:1.; double N_e = (1.-tau)*N_0;
	    double temp    = esyn/cst::mec2;
	    if(temp>0.0) temp=1./temp; 
	    double y = (*x1)*temp; //adim
	    double f1= pow(y,0.297)*exp(-y); //adim 
	    
	    if (gi<gamma_min ) 
	      {
		N_e *= gi/gamma_min*pow(gamma_min,-cst::p); //adim;
	      } 
	    else if (gi>gamma_max) 
	      {
		N_e *= 0.0;
	      }
	    else 
	      {
		N_e *= pow(gi,-cst::p); //adim;
	      }	      
	    
	    
	    if(ComovingTime<tcross)
	      {
		N_e *= tsyn/tcross*(1.-exp(-ComovingTime/tsyn));
	      }
	    else
	      {
		N_e *= tsyn/tcross*(1.-exp(-tcross/tsyn))*exp((tcross-ComovingTime)/tsyn);
	      }
	    
	    (*its)  = f1*N_e*1.104; //adim
	    
	    its++;
	    x1++;
	  }	      
	//OBSERVED SPECTRUM
	// P is a power -> P_obs = P_em
	// Esyn         -> LT * Esyn
	// Psyn/Esyn    -> Psyn/Esyn/LT [MeV/MeV/s]
	// Psyn/Esyn/DE -> Psyn/Esyn/LT^2 [1/MeV/s]
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
  }
}

