#include "GRBSynchrotron.h"
#include "GRBConstants.h"
#include "SpectObj.h"


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
  if (time<=0.) 
    {
      m_spectrumObj *= 0.0; 
      return;
    }
  const int type=0; // Power Law
  //  const int type=1; // Complete integration;
  double ComovingTime; 
  const double EnergyTransformation = 
    m_spectrumObj.getEnergyTransformation(GAMMAF,angle); // in eV
  //m_spectrumObj.obs2em(GAMMAF,angle);
  
  //////////////////////////////////////////////////
  // Gamma e+-
  //////////////////////////////////////////////////
  
  const double Psyn_0 = 4./3.*Umag(B)*cst::erg2MeV*cst::c*cst::st; //MeV/s
  const double esyn_0 = (3.*cst::pi/8.)*(B/cst::BQ)*cst::mec2; //MeV
  const double Psyn_esyn = Psyn_0/esyn_0; //1/sec
  const double tau_0 = (cst::flagIC == 1.) ? cst::c*cst::st*N0/Vol 
    : cst::flagIC; //1/s
  const double tcross  = dr/cst::c;
  
  //////////////////////////////////////////////////
  /*
    std::vector<double>::iterator x1 = x.begin();
    std::vector<double>::iterator its = fsyn.begin();
  */
  double   dgamma     = pow((gamma_max/gamma_min),1./cst::enstep);
  int it,j;
  for(it = 0; it < m_spectrumObj.size();) 
    {
      double value=0.0;
      // Observed Energy, in eV
      double obs_energy = m_spectrumObj.getEnergy(it);
      ComovingTime = comovingTime(time,GAMMAF,obs_energy,
				  distance_to_source);
      //comoving time <0 means that dispersive effect delayed 
      //to much the photon at this energy bin
      if(ComovingTime<=0){it++; continue;}
      
      // Observed Energy, in mec2 units
      obs_energy *= (1.e-6/cst::mec2);
      // Emitted  Energy, in mec2 units
      double em_energy  = obs_energy/EnergyTransformation;
      
      if(type == 0) //Power Law approx...
	{
	  double gc = cst::mec2/(Psyn_0*ComovingTime); 
	  gc = (gc <= 1. ? 1.+1.e-6 : gc);
	  double em = 4./3.*esyn_0*(pow(gamma_min,2.)-1.)/cst::mec2;
	  // The SYN energy that conrespond to gc 
	  double ec = 4./3.*esyn_0*(pow(gc,2.)-1.)/cst::mec2;
	  double gi      = 2.*gamma_min; 
	  double gi2     = pow(gi,2.);
	  double tsyn    = (gi*cst::mec2)/(Psyn_0*(gi2-1.));  
	  double tau = tau_0;
	  tau *= (tsyn<tcross)? tsyn: tcross;
	  if(tau>1.) tau = 1.;
	  double N_e = (1.-tau)*electronNumber(gi, gamma_min, gamma_max,
					       dr, ComovingTime, tsyn, N0);
	  
	  //	  if(*x1<gc) //ekn == gc
	  value = processFlux(em_energy,ec,em) //adim
	    *N_e*gi*(dgamma-1.); 
	}
      else // Complete integration
	{
	  // Here the integral beween the population of electron and 
	  // the syncrotron emission from each of those electrons.
	  for (j = 0; j < cst::enstep;)
	    {
	      double gi      = gamma_min*pow(dgamma,j);
	      double gi2     = pow(gi,2.);
	      double esyn    = esyn_0*(gi2-1.)/cst::mec2; // In mec2 units
	      // an electronb with gi, emits its energy e0 = (gi*cst::mec2) with a power
	      // Psyn = Psyn_0*(gi2-1.), in a time tsyn = e0/Psyn
	      double tsyn    = (gi*cst::mec2)/(Psyn_0*(gi2-1.)); 
	      double tau = tau_0; //1/sec
	      tau *= (tsyn<tcross)? tsyn : tcross;
	      if(tau>1.) tau = 1.;

	      /*
		double N_e = electronNumber(gi, gamma_min, gamma_max,
		N0);
	      */
	      //Test the number of electron:
	      /*
		value += N_e*gamma_min*(pow(dgamma,j+1)-pow(dgamma,j));
	      */
	      //Test the total electron energy: 0.895193
	      //	      value += cst::mec2*gi*N_e*gamma_min*(pow(dgamma,j+1)-pow(dgamma,j));
	      //Test the temporal evolution:
	      /*
		value += TemporalEvolution(gi,
		gamma_min,
		dr,
		ComovingTime, 
		tsyn); // *gamma_min*(pow(dgamma,j+1)-pow(dgamma,j));
	      */
	      double dg = gamma_min*(pow(dgamma,j+1)-pow(dgamma,j)); //adim 
	      double N_e = (1.-tau)
		*electronNumber(gi, gamma_min, gamma_max,
				dr, ComovingTime, tsyn, N0);
	      value+=N_e*SynchrotronFunction(esyn,em_energy)*dg;
	      
	      
	      j+=1; // This is for speed up the simulation
	    }
	}
      //      cout<<(value/cst::erg2MeV)/(Umag(B)*Vol)<<" "<<Umag(B)*Vol<<endl;
      //cout<<value<<" "<<endl;
      m_spectrumObj.SetSpectrum(it,value*1.e+6/m_spectrumObj.getBin(it) );
      // 1/MeV
      it++;
    } 
  
  m_spectrumObj *= (Psyn_esyn); // 1/MeV/s
  m_spectrumObj.em2obs(GAMMAF,angle);  // is in ph/s/MeV
  //  return m_spectrumObj; 
} 

double GRBSynchrotron::SynchrotronFunction(double esyn, double energy)
{
  // note that energy and esyn are in unit of mec2
  if(esyn <= 0.0) return 0.0;
  double y       = energy/esyn; //adim
  double f1      = pow(y,0.297)*exp(-y); //adim
  return f1;
}

