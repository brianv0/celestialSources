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
  if (time<=0.) {
    m_spectrumObj *= 0.0; 
    /*
      cout<<"   SYNCRO      "<<endl;
      cout<<"   Rise Time  = "<<(dr/cst::c)/GAMMAF<<endl;
      cout<<"   Decay Time = "<<gamma_min*cst::mec2/
      (4./3.* Umag(B)*cst::erg2MeV*cst::c*cst::st*
      (pow(gamma_min,2.)-1.))/GAMMAF<<endl;
      cout<<" Gamma_min    = "<<gamma_min<<endl;
    */
    return;}
  
  double ComovingTime; 
  const double EnergyTransformation = 
    1.e+6 * cst::mec2 *
    ( GAMMAF + sqrt(GAMMAF*GAMMAF+1.)*cos(angle) );
  
  //////////////////////////////////////////////////
  // Gamma e+-
  //////////////////////////////////////////////////
  
  std::vector<double>  x  = 
    m_spectrumObj.getEnergyVector(1./(EnergyTransformation));
  std::vector<double>  dx = 
    m_spectrumObj.getBinVector(1./(EnergyTransformation));
  std::vector<double> fsyn(cst::enstep+1, 0.0);
  
  const double Psyn_0 = 4./3.*Umag(B)*cst::erg2MeV*cst::c*cst::st; //MeV/s
  const double esyn_0 = (3.*cst::pi/8.)*(B/cst::BQ)*cst::mec2; //MeV
  const double Psyn_esyn = Psyn_0/esyn_0; //1/sec
  const double tau_0 = (cst::flagIC == 1.) ? cst::c*cst::st*N0/Vol 
    : cst::flagIC; //1/s
  const double tcross  = dr/cst::c;
  
  //////////////////////////////////////////////////
  
  std::vector<double>::iterator x1 = x.begin();
  std::vector<double>::iterator its = fsyn.begin();
  double   dgamma     = pow((gamma_max/gamma_min),1./cst::enstep);
  

  
  while(x1!=x.end()) 
    {
      ComovingTime = comovingTime(time,GAMMAF,EnergyTransformation*(*x1),
				  distance_to_source);
      //comoving time <0 means that dispersive effect delayed 
      //to much the photon at this energy bin
      if(ComovingTime<=0){its++; x1++; continue;}
      
      // Here the integral beween the population of electron and 
      // the syncrotron emission from each of those electrons.
      for (int j = 0; j < cst::enstep;)
	{
	  double gi      = gamma_min*pow(dgamma,j);
	  double gi2     = pow(gi,2.);
	  double esyn    = esyn_0*(gi2-1.)/cst::mec2; // In mec2 units
	  
	  double tsyn    = (gi*cst::mec2)/(Psyn_0*(gi2-1.)); 
	  double tau = tau_0;
	  tau *= (tsyn<tcross)? tsyn: tcross;
	  if(tau>1.) tau = 1.;
	  double N_e     = (1.-tau)*electronNumber(gi, gamma_min, gamma_max,
						   dr, ComovingTime, tsyn, N0);
	  (*its) += SynchrotronFunction(esyn,(*x1))*N_e*(pow(dgamma,j+2)-pow(dgamma,j)); //adim
	  //cout<<" SYN = "<<(*its)<<endl;
	  j+=2; // This is for speed up the simulation
	}
      its++;
      x1++;
    } 
  //OBSERVED SPECTRUM
  // LT = Lorentz Transformation
  // P  is a power -> P_obs = P_em
  // Esyn          -> LT * Esyn
  // Psyn/Esyn     -> Psyn/Esyn/LT [MeV/MeV/s]
  // Psyn/Esyn/DE  -> Psyn/Esyn/LT^2 [1/MeV/s]
  
  m_spectrumObj.clear();
  x1  =    x.begin();
  its = fsyn.begin();
  std::vector<double>::iterator dx1 =   dx.begin();
  
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

double GRBSynchrotron::SynchrotronFunction(double esyn, double energy)
{
  // note that energy and esyn are in unit of mec2
  esyn = (esyn>0.0) ? 1./esyn : 0.0; 
  double y       = energy*esyn; //adim
  double f1      = pow(y,0.297)*exp(-y); //adim
  return f1;
}

