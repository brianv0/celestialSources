#include <iostream>
#include <algorithm>

#include "GRBengine.h"
#include "GRBShock.h"
#include "GRBShell.h"

#define DEBUG 0

using namespace cst;

double Beta(double gamma)
{
  return  sqrt(gamma*gamma-1.)/gamma;
}

GRBengine::GRBengine(Parameters *params)
  : m_params(params)
{
  m_dir = m_params->GetGalDir();
  
  std::cout<<" Create new GRB : "<<std::endl;
  std::cout<<" Seed           = "<<m_params->GetGRBNumber()<<std::endl;
  std::cout<<" Position    : l="<<m_dir.first<<", b = "<<m_dir.second<<std::endl;
  
}

std::vector<GRBShock*> GRBengine::CreateShocksVector()
{
  //////////////////////////////////////////////////
  const   int    Nshell            = m_params->GetNshell();
  const   double Etot              = m_params->GetEtot();
  const   double InitialSeparation = m_params->GetInitialSeparation();
  const   double InitialThickness  = m_params->GetInitialThickness();
  const   double Gmin              = m_params->GetGammaMin();
  const   double Gmax              = m_params->GetGammaMax();
  std::vector<GRBShell*> theShells;
  std::vector<GRBShock*> theShocks;
  
  double tobs   = 0.0;
  double tGRB   = 0.0;
  
  for(int i = 0; i < Nshell ; i++ )
    {
      double g;
      if(Nshell ==2 && i==1 )
	g += m_params->rnd->Uniform(Gmax-g+1,Gmax); 
      else 
	g = m_params->rnd->Uniform(Gmin,Gmax); // lorentz factor;
      double r  = InitialSeparation + (Nshell-(i+1))*(InitialSeparation+InitialThickness); // initial radius
      //      double dr = r/(g*g); 
      GRBShell *ashell;
      ashell = new GRBShell(g,r,InitialThickness,Etot);
      theShells.push_back(ashell);
    }
  
  //  int N      = theShells.size();
  
  bool Continue = true;

  while (theShells.size()>1)
    {
      double ShockTime     = 1e50;
      int FirstShell       =  -1;
      GRBShell *sh1,*sh2;
      for(int i = 0; i < theShells.size()-1 ; i++ )
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

      if(FirstShell>-1)
	{  
	  for(int i = 0; i < theShells.size(); i++ )
	    {  
	      theShells[i]->Evolve(ShockTime);
	    }
	  
	  tGRB += ShockTime;
	  sh1  = theShells[FirstShell];
	  sh2  = theShells[FirstShell+1];
	  
	  tobs = tGRB - sh2->GetRadius()/c;
	  if(DEBUG) std::cout<<" Shock: GRB time = "<<tGRB<<" Obs time ="<<tobs<<" radius = "<<sh2->GetRadius()<<std::endl;
	  
	  GRBShock *ashock = new GRBShock(sh1,sh2,tobs);
	  theShocks.push_back(ashock);	 
	  theShells[FirstShell] = ashock->MergedShell();
	  std::vector<GRBShell*>::iterator pos = theShells.begin();
	  for (int i = 0;i<FirstShell+1;i++) pos++;
	  theShells.erase(pos);
	}
      else
	if(theShocks.size()==0) 
	  {
	    std::cout<<" RECOMPUTING ShocksVector "<<std::endl;
	    for(int i =0; i<(int) theShells.size();i++)
	      delete theShells[i];
	    return CreateShocksVector();
	  }
	else
	  break;
    }
  //////////////////////////////////////////////////
  
  for(int i =0; i<(int) theShells.size();i++)
    delete theShells[i];
  
  std::vector<GRBShock*>::iterator pos = theShocks.begin();
  
  while(pos!=theShocks.end())
    {
      if((*pos)->GetEfficiency()<1.0e-2)
	theShocks.erase(pos);
      else
	{
	  pos++;
	}
    }

  if(theShocks.size()==0) 
    {
      return CreateShocksVector();
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
      for(int i = 0; i< (int) theShocks.size(); i++)
	{   
	  std::cout<<"---- Shock N = "<<i
		   <<" at tobs = "<<theShocks[i]->GetTime()
		   <<" Efficiency: "<<theShocks[i]->GetEfficiency()
		   <<std::endl;
	}
    }
  return theShocks;
}

//////////////////////////////////////////////////

/*
  double GRBengine::getDistance()
  {
  double m_redshift = cst::red;
  double qo=(1.0+3.0*cst::wzel)/2.0;
  return ((cst::c/(cst::Hubble*1.0e+5)/pow(qo,2.0))*
  (m_redshift*qo+(qo-1.0)*(-1.0+sqrt(2.0*qo*m_redshift+1.0)))*cst::mpc2cm);
  }
*/


