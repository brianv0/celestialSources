#include <iostream>
#include <algorithm>

#include "GRBengine.h"
#include "GRBShock.h"
#include "GRBShell.h"

#define DEBUG 0

using namespace cst;

GRBengine::GRBengine(Parameters *params)
  : m_params(params)
{
  m_dir = m_params->GetGalDir();
}

std::vector<GRBShock*> GRBengine::CreateShocksVector()
{
  //////////////////////////////////////////////////
  double tobs=0.0;
  
  const   double Gmin              = m_params->GetGammaMin();
  const   double Gmax              = m_params->GetGammaMax();
  const   double InitialSeparation = m_params->GetInitialSeparation();
  const   double InitialThickness  = m_params->GetInitialThickness();
  const   double Etot              = m_params->GetEtot();
  const   double BurstDuration     = m_params->GetDuration();  
  

  GRBShock *ashock;
  std::vector<GRBShock*> theShocks;
  double tshock,rshock,betaSlow;
  int i=0;
  double Gmean=0.0;
  
  while(tobs < BurstDuration)
    {
      double G1 = m_params->rnd->Uniform(Gmin,Gmax); // lorentz factor;
      double G2 = m_params->rnd->Uniform(Gmin,Gmax); // lorentz factor;
      double GammaSlow = TMath::Min(G1,G2);
      double GammaFast = TMath::Max(G1,G2);

      double M1 = Etot/(GammaSlow*c2);
      double M2 = Etot/(GammaFast*c2);
      double mf  = M1 + M2; //g
      double M   = sqrt(M1*M1+M2*M2+2.0*M1*M2*(GammaSlow*GammaFast-sqrt(GammaSlow*GammaSlow-1.0)*sqrt(GammaFast*GammaFast-1.0)))-(mf);
      
      double GammaMean = (M1*GammaSlow + M2*GammaFast)/(mf+M);//0.5*(GammaFast+GammaSlow);
      Gmean+=GammaMean;
      
      rshock    = pow(GammaMean,2.) * InitialSeparation;
      betaSlow  = sqrt(1.0- 1.0/pow(GammaSlow,2.0));
      tshock    = 2.0*i*m_params->rnd->Uniform(0,1)*InitialSeparation/c + rshock/(c * betaSlow);
      
      GRBShell *FastShell = new GRBShell(GammaFast,rshock,InitialThickness,Etot);
      GRBShell *SlowShell = new GRBShell(GammaSlow,rshock,InitialThickness,Etot);
      
      ashock    = new GRBShock(SlowShell,FastShell,tshock,2.5);//m_params->rnd->Uniform(2.0,3.0));
      theShocks.push_back(ashock);

      delete FastShell;
      delete SlowShell;

      std::sort(theShocks.begin(), theShocks.end(), ShockCmp());
      tobs=theShocks[0]->GetTime()+theShocks.back()->GetTime()+2.0*theShocks.back()->GetDuration();
      
      
      if(DEBUG) 
	{
	  
	  double d7    = InitialThickness*1.0e-7;
	  double r10   = rshock*1e-10;
	  double g100  = GammaMean/100.0;
	  double rsh14 = rshock * 1.0e-14; 
	  double E52   = Etot * 1e-52;
	  double n1    = (Etot * erg2meV)/(GammaMean * mpc2 * 4.0 *pi *rshock * rshock * InitialThickness * GammaMean);
	  
	  double Ep  = 363.742   * sqrt(E52/d7 * ab*3.0)*pow(ae*3.,2.)/rsh14;
	  double B   = 5.13809e6 * sqrt(E52/d7 * ab*3.0)              /(rsh14*g100); 
	  
	  std::cout<<" Shock: GRB time = "<<
	    tshock<<" Obs time ="<<ashock->GetTime()<<" radius = "<<ashock->GetRadius()<<" "<<rshock<<" efficiency : "<<ashock->GetEfficiency()<<" Ep = "<<Ep<<std::endl;
	  std::cout<<" Shocks size = "<<theShocks.size()<<std::endl;	    
	}
      i++;
    }
  std::sort(theShocks.begin(), theShocks.end(), ShockCmp());
  double T0 = theShocks[0]->GetTime();
  for(int i = 1; i<= (int) theShocks.size(); i++)
    {
      theShocks[theShocks.size()-i]->
	SetTime(theShocks[theShocks.size()-i]->GetTime() - T0);
    }

  if(DEBUG)
    {
      int i=0;
      for(i = 0; i< (int) theShocks.size(); i++)
	{   
	  theShocks[i]->Print();
	}
      std::cout<<i<<" shocks computed..."<<std::endl;
    }
  
  return theShocks;
}

//////////////////////////////////////////////////

double GRBengine::GetDistance()
{
  static const double z = 1.0;
  static const double Hubble = 75.0 * (3.24e-20); // 1/s
  double W = 1.0;
  double Wm = 0.3;
  double Wl = W-Wm;
  double qo = (W/2. - Wl);
  double dL = 2.0*cst::c/(Hubble*sqrt(Wm))*(1.-pow(1.+z,-1./2.));
  return dL;
}


