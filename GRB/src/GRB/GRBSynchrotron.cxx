#include "GRBSynchrotron.h"
#include "GRBConstants.h"

GRBSynchrotron::GRBSynchrotron()
  : RadiationProcess()
{}

GRBSynchrotron::GRBSynchrotron(SpectObj spectrumObj)
  : RadiationProcess(spectrumObj)
{}

void GRBSynchrotron::load(const GRBShock *Shock, 
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

void GRBSynchrotron::load(const double time,
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
  if (time<=0.) {m_spectrumObj *= 0.0; return;}
  
  int type = 0; // 1-> integration
                // 0-> Power Law
  double ComovingTime; 
  const double EnergyTransformation = 
    1.e+6 * cst::mec2 *
    ( GAMMAF + sqrt(GAMMAF*GAMMAF+1.)*cos(angle) );
  
  //////////////////////////////////////////////////
  // Gamma e+-
  //////////////////////////////////////////////////
  
  std::vector<double>  x  = m_spectrumObj.getEnergyVector(1./(EnergyTransformation));
  std::vector<double>  dx = m_spectrumObj.getBinVector(1./(EnergyTransformation));
  std::vector<double> fsyn(cst::enstep+1, 0.0);
  
  const double Psyn_0 = 4./3. * Umag(B) * cst::erg2MeV * cst::c * cst::st; //MeV/s      
  const double esyn_0 = (3.*cst::pi/8.)*(B/cst::BQ)*cst::mec2; //MeV
  const double Psyn_esyn = Psyn_0/esyn_0; //1/sec
  const double tau_0 = (cst::flagIC == 1.) ? cst::c*cst::st*N0/Vol : cst::flagIC; //1/s
  
  if (type==0) // POWER LAW
    {
      std::vector<double>::iterator x1 = x.begin();
      std::vector<double>::iterator its = fsyn.begin();
      while(x1!=x.end()) 
	{
	  ComovingTime = comovingTime(time,GAMMAF,(*x1),distance_to_source);
	  //comoving time <0 means that dispersive effect delayed 
	  //to much the photon at this energy bin
	  if(ComovingTime<=0){its++; x1++; continue;}

	  double gc = cst::mec2/(Psyn_0*ComovingTime);
	  gc = (gc <= 1. ? 1.+1.e-6 : gc);
	  
	  //em and ec are specific:
	  double em = esyn_0*(pow(gamma_min,2.)-1.)/cst::mec2;
	  double ec = esyn_0*(pow(gc,2.)-1.)/cst::mec2;
	  double gi = sqrt(1. + (*x1)*cst::mec2/esyn_0);
	  //cout<<em*1.0e-6*EnergyTransformation<<" "<<ec*1.0e-6*EnergyTransformation<<endl;
	  double gi2     = pow(gi,2.);
	  double Psyn    = Psyn_0*(gi2-1.); //MeV/s
	  double tsyn    = (gi*cst::mec2)/Psyn;
	  double tsyn0    = (cst::mec2)/(gamma_min*Psyn_0);
	  double tcross  = dr/(cst::c);
	  double tau = tau_0;
	  tau *= (tsyn<tcross)? tsyn: tcross;
	  if(tau>1.) tau = 1.;
	  
	  double N_e = (1.- tau)*N0;	  
	  
	  if(ComovingTime <= tcross)
	    N_e *= tsyn0/tcross*(1.-exp(-ComovingTime/tsyn0));
	  else
	    N_e *= tsyn0/tcross*(1.-exp(-tcross/tsyn0))*exp((tcross-ComovingTime)/tsyn0);
	  
	  (*its) = processFlux(*x1,ec,em);
	  (*its++) *= N_e; // adim
	  (x1++);
	}
    }
  else
    {
      //////////////////////////////////////////////////
      // Integration Needs to  be fixed !!            //
      // The integration was corrected :(
      ////////////////////////////////////////////////// 	  
      
      std::vector<double>::iterator x1 = x.begin();
      std::vector<double>::iterator its = fsyn.begin();
      double   dgamma     = pow((gamma_max/gamma_min),1./cst::enstep);
      
      while(x1!=x.end()) 
	{
	  ComovingTime = comovingTime(time,GAMMAF,(*x1),distance_to_source);
	  //comoving time <0 means that dispersive effect delayed 
	  //to much the photon at this energy bin
	  if(ComovingTime<=0){its++; x1++; continue;}
	  
	  for (int j = 0; j < cst::enstep; j++)
	    {
	      double gi      = gamma_min*pow(dgamma,j);
	      double gi2     = pow(gi,2.);
	      double esyn    = esyn_0*(gi2-1.); //MeV
	      double temp    = esyn/cst::mec2;
	      temp = (temp>0.0) ? 1./temp : 0.0; 
	      double y = (*x1)*temp; //adim
	      double f1= pow(y,0.297)*exp(-y); //adim 
	      double N_e = electronNumber(gi, gamma_min, gamma_max,
					  dr, ComovingTime, Psyn_0, N0);
	      (*its) += f1*N_e*gamma_min*(pow(dgamma,j+1)-pow(dgamma,j));//dg[j]; //adim
	    }
	  
	  its++;
	  x1++;
	}
      
      
      //OBSERVED SPECTRUM
      // P  is a power -> P_obs = P_em
      // Es yn         -> LT * Esyn
      // Psy n/Esyn    -> Psyn/Esyn/LT [MeV/MeV/s]
      // Psyn /Esyn/DE -> Psyn/Esyn/LT^2 [1/MeV/s]
    } 
  //// //////////////////////////////////////////////
  // COMMON 
  //////////////////////////////////////////////////
  m_spectrumObj.clear();
  std::vector<double>::iterator x1  =    x.begin();
  std::vector<double>::iterator dx1 =   dx.begin();
  std::vector<double>::iterator its = fsyn.begin();
  while(dx1!= dx.end()) 
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
  m_spectrumObj.SetSpectrum(x,fsyn);  // is in ph/s/MeV
  return; 
} 


