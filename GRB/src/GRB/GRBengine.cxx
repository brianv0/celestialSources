#include "GRBengine.h"
#include "GRBShock.h"
#include "GRBShell.h"

using namespace cst;

double Beta(double gamma)
{
  return  sqrt(gamma*gamma-1.)/gamma;
}

GRBengine::GRBengine(Parameters *params)
  : m_params(params)
{
  
  double l = m_params->rnd->Uniform(180.0,-180.0);
  double b = m_params->rnd->Uniform(90.0,-90.0);
  
  m_dir = std::make_pair(l,b);
  
  std::cout<<" Create new GRB : "<<std::endl;
  std::cout<<" Seed           = "<<m_params->GetGRBNumber()<<std::endl;
  std::cout<<" Position    : l="<<l<<", b = "<<b<<std::endl;
  
}

std::vector<GRBShock*> GRBengine::CreateShocksVector(const int Nshell, 
						     const double L, //initial separation 
						     const double l, //initial tickness
						     const double etot)
{
  std::vector<GRBShell*> theShells;
  std::vector<GRBShock*> theShocks;
  
  double tobs   = 0.0;
  double tGRB   = 0.0;
  
  for(int i = 0; i < Nshell ; i++ )
    {
      double g  = m_params->rnd->Uniform(gmin,gmax); // lorentz factor;
      double r  = L + (Nshell-(i+1))*(L+l); // initial radius
      //      double dr = r/(g*g); 
      GRBShell *ashell = new GRBShell(g,r,l,etot/Nshell);
      theShells.push_back(ashell);
    }
  
  int N      = theShells.size();
  int Nshock = 0;
  while (N > 1)
    {
      double ShockTime     = 1e50;
      int FirstShell       =  -1;
      GRBShell *sh1,*sh2;
      for(int i = 0; i < N-1 ; i++ )
	{
	  sh1 = theShells[i];
	  sh2 = theShells[i+1];
	  double b1 = sh1->GetBeta();
	  double b2 = sh2->GetBeta();
	  double D  = sh1->GetRadius() - (sh2->GetRadius()+ sh2->GetThickness());
	  
	  if (b1 < b2)
	    { 
	      if(D/(c*(b2-b1)) <= ShockTime)
		{
		  ShockTime  = D/(c*(b2-b1)); 
		  FirstShell = i;
		}
	    }
	}
      if(FirstShell==-1) break;
      
      for(int i = 0; i < N; i++ )
	{  
	  theShells[i]->Evolve(ShockTime);
	}
      
      tGRB += ShockTime;
      sh1  = theShells[FirstShell];
      sh2  = theShells[FirstShell+1];
      
      tobs = tGRB - sh2->GetRadius()/c;
      //std::cout<<" Shock: GRB time = "<<tGRB<<" Obs time ="<<tobs<<" radius = "<<sh2->GetRadius()<<std::endl;
      
      GRBShock *ashock = new GRBShock(sh1,sh2,tobs);
      theShocks.push_back(ashock);	 
      GRBShell *Ms = ashock->MergedShell();
      theShells[FirstShell] = Ms;
      std::vector<GRBShell*>::iterator pos = theShells.begin();
      for (int i = 0;i<FirstShell+1;i++) pos++;
      theShells.erase(pos);
      Nshock++;
      N = theShells.size();
      //std::cout<<"Nshock = "<<Nshock<<" Nshell "<<N<<std::endl;
    }
  
  std::sort(theShocks.begin(), theShocks.end(), ShockCmp());
  double T0 = theShocks[0]->GetTime();
  
  for(int i = 1; i<= (int) theShocks.size(); i++)
    {
      theShocks[theShocks.size()-i]->
	SetTime(theShocks[theShocks.size()-i]->GetTime() - T0);
    }
  
  for(int i = 0; i< (int) theShocks.size(); i++)
    {   
      std::cout<<"---- Shock N = "<<i
	       <<" at tobs = "<<theShocks[i]->GetTime()
	       <<std::endl;
      //      GRBShell *Ms = theShocks[i]->MergedShell();
    }
  return theShocks;
}

//////////////////////////////////////////////////
double GRBengine::GenerateGamma(double gammamin,double gammamax)
{
  // cout<<1.+fabs(m_params->rnd->Gaus(gammamax-gammamin,gammamax-gammamin))<<endl;
  // return 1.+fabs(m_params->rnd->Gaus(gammamax-gammamin,gammamax-gammamin));
  // return 1.+m_params->rnd->Exp(gammamax-gammamin/2);//,gammamax-gammamin));
  return m_params->rnd->Uniform(gammamin,gammamax);
}

/*
  double GRBengine::getDistance()
  {
  double m_redshift = cst::red;
  double qo=(1.0+3.0*cst::wzel)/2.0;
  return ((cst::c/(cst::Hubble*1.0e+5)/pow(qo,2.0))*
  (m_redshift*qo+(qo-1.0)*(-1.0+sqrt(2.0*qo*m_redshift+1.0)))*cst::mpc2cm);
  }
*/


