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
  double tv     = (InitialSeparation+InitialThickness)/c;
  //////////////////////////////////////////////////
  // 1) Create the shells...
  // )    )       )    )    )
  // N   ...      2    1    0
  double g,b,r;
  double temission=0;
  GRBShell *ashell;
  GRBShock *ashock;
  
  while(theShells.size()< Nshell)
    {
      //      g = m_params->rnd->Uniform(Gmin,Gmax); // lorentz factor;
      //g = m_params->rnd->Uniform(Gmin,Gmax); // lorentz factor;
      //g = TMath::Max(1.0,m_params->rnd->Exp((Gmax+Gmin)/2.0)); // lorentz factor;
      g = TMath::Max(1.0,m_params->rnd->Gaus((Gmax+Gmin)/2.0,(Gmax-Gmin)/4.0)); // lorentz factor;
      b =  sqrt(1.0- 1.0/pow(g,2.0));
      //      if(m_params->rnd->Uniform() < m_params->GetTau())
      //	{      
      //      r  = InitialSeparation + temission * b * c;
      ashell = new GRBShell(g,0,InitialThickness,Etot,temission);
      theShells.push_back(ashell);//insert(theShells.begin(),ashell);
      //      theShells.insert(theShells.begin(),ashell);
      //}
      temission+=tv;
      if(DEBUG) std::cout<<" new Shell, shell size = "<<theShells.size()<<std::endl;
      
    }
  
  if(DEBUG)
    std::cout<<" Shells generated: shell size = "<<theShells.size()<<" N shell = "<<Nshell<<std::endl;
  
  double b1, b2,D;
  double r1, r2, dr2;
  double t1,t2,tsh;
  while (theShells.size()>1)
    {
      double ShockTime     = 1e50;
      int FirstShell       =  -1;
      GRBShell *sh1,*sh2;
      
      for(int i = 0; i < theShells.size()-1 ; i++ )
	{
	  sh1 = theShells[i];
	  sh2 = theShells[i+1];
	  b1  = sh1->GetBeta();
	  b2  = sh2->GetBeta();
	  t1  = sh1->GetTime();
	  t2  = sh2->GetTime();
	  if(t2 > t1 && b2 > b1)
	    {
	      dr2 = sh2->GetThickness();
	      r2  = sh2->GetRadius(t2);
	      r1  = sh1->GetRadius(t2);
	      D   = r1 + b2*(t2-t1)*c -(r2+dr2);
	      tsh = t2 + D/(c*(b2-b1));
	      
	      if(tsh <= ShockTime)
		{
		  ShockTime  = tsh;
		  FirstShell = i;
		}
	    }
	}
      if(FirstShell>-1)
	{  
	  tGRB = ShockTime;
	  sh1  = theShells[FirstShell];
	  sh2  = theShells[FirstShell+1];
	  ashock = new GRBShock(sh1,sh2,tGRB,2.5);//m_params->rnd->Uniform(2.0,3.0));

	  theShocks.push_back(ashock);
	  theShells[FirstShell] = ashock->MergedShell();

	  std::vector<GRBShell*>::iterator shelliter = theShells.begin();
	  for (int i = 0;i<FirstShell+1;i++) shelliter++;
	  theShells.erase(shelliter);
	  
	  if(DEBUG) 
	    {
	      std::cout<<" Shock: GRB time = "<<
		tGRB<<" Obs time ="<<ashock->GetTime()<<" radius = "<<ashock->GetRadius()<<" efficiency : "<<ashock->GetEfficiency()<<std::endl;
	      std::cout<<" Shells size = "<<theShells.size()<<" Shocks size = "<<theShocks.size()<<std::endl;	    
	    }
	}
      else
	if(theShocks.size()==0) 
	  {
	    if(DEBUG) 
	      std::cout<<" theShocks.size() == 0 => RECOMPUTING ShocksVector "<<std::endl;
	    theShells.erase(theShells.begin(),theShells.end());
	    return CreateShocksVector();
	  }
      	else
	  {
	    break;
	  }
    }

  //////////////////////////////////////////////////
  theShells.erase(theShells.begin(),theShells.end());
  
  std::vector<GRBShock*>::iterator pos = theShocks.begin();
  
  while(pos!=theShocks.end())
    {
      if((*pos)->GetEfficiency() < 1.0e-3 || (*pos)->GetRadius() > 1e17)
	theShocks.erase(pos);
      else
	{
	  pos++;
	}
    }
  
  if(theShocks.size()==0) 
    {
      theShells.erase(theShells.begin(),theShells.end());
      
      if(DEBUG)
	std::cout<<" theShocks.size()==0  =>   RECOMPUTING THE SHOCKS VETOR = "<<std::endl;
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
		   <<" radius = "<<theShocks[i]->GetRadius()
		   <<std::endl;
	}
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
  double dL = (cst::c/Hubble) * z + pow(z,2.0)/2.*(1.0-qo);
  return dL;
}


