#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

#include "GRBConstants.h"

#define DEBUG  0

Parameters::Parameters()
{
  rnd = new TRandom();
  SetGRBNumber(UInt_t(rnd->GetSeed()));
  m_Type=0;
  SetGBMOutput(false);
}

void Parameters::SetGRBNumber(UInt_t GRBnumber)
{
  m_GRBnumber = GRBnumber;
  rnd->SetSeed(m_GRBnumber);
  double tmp = rnd->Uniform();
}

//////////////////////////////////////////////////
double Parameters::GetBATSEFluence()
{
  using std::pow;
  if (m_Type==1)
    return pow(10.0,(double)rnd->Gaus(-6.3,0.57)); //erg/cm^2 (Short Bursts)
  return pow(10.0,(double)rnd->Gaus(-5.4,0.62)); //erg/cm^2 (Long Burst)
}

//////////////////////////////////////////////////
void Parameters::SetGalDir(double l, double b)
{
  double ll,bb;
  
  ll = (l<=180.0 && l>=-180.0) ? l : rnd->Uniform(-180.0,180.0);
  bb = (b<=90.0 && b>=-90.0)   ? b : rnd->Uniform(-90.0,90.0);
  
  m_GalDir=std::make_pair(ll,bb);
}

void Parameters::SetNshell(int nshell)
{
  if(nshell>1 && nshell<4)
    m_Type=1;
  else if(nshell>4)
    m_Type=2;
  
  if (nshell<1)
    if (m_Type==1)    
      m_Nshell = int(2+rnd->Integer(2));
    else
      m_Nshell = int(10+rnd->Integer(70));
  else
    m_Nshell = nshell;
}

void Parameters::SetFluence(double fluence)
{
  if(fluence == 0)  
    m_Fluence = GetBATSEFluence();
  else
    {
      m_Fluence = fluence;
      double tmp = rnd->Uniform();//THB: this needed parentheses
    }
}

void Parameters::SetEtot(double etot)
{
  m_Etot = (etot>0.0) ? etot : pow(10.0,rnd->Gaus(51.,0.5)); //erg
}

/*
  void Parameters::SetEpeak(double epeak)
{
  if(epeak<=0)
    {
      if (m_Type==1)
	m_Ep = pow(10.0,(double)rnd->Gaus(2.0,0.5)); //erg/cm^2 (Short Bursts)
      else
	m_Ep = pow(10.0,(double)rnd->Gaus(2.5,0.5)); //erg/cm^2 (Long Burst)
    }
  else 
    m_Ep =  epeak;
}
*/

void Parameters::SetInitialSeparation(double initialSeparation)
{
  m_InitialSeparation = (initialSeparation>0.0) ? initialSeparation : pow(10.,rnd->Uniform(7.0,10.0));
}

void Parameters::SetInitialThickness(double initialThickness)
{
  m_InitialThickness = (initialThickness>0.0) ? initialThickness : pow(10.,rnd->Uniform(7.0,10.0));
}
  
void Parameters::SetGammaMin(double gmin)
{
  m_Gmin = (gmin>=1.0) ? gmin : 100.0;
}

void Parameters::SetGammaMax(double gmax)
{
  m_Gmax = (gmax>=1.0 && gmax>=m_Gmin) ? gmax : 10.0*m_Gmin;
}

void Parameters::SetInverseCompton(double ic)
{
  m_InverseCompton = (ic>=0.0) ? ic : rnd->Uniform(0.,10.);
}

//..................................................

void Parameters::ReadParametersFromFile(std::string paramFile, int NGRB)
{
  
  std::ifstream f1(paramFile.c_str());
  if (!f1.is_open()) 
    {
      std::cout<<"GRBConstants: Error Opening paramFile!! \n";
      exit(1);
    }
  int    nshell;
  double fluence,etot,r0,dr0;
  double gmin,gmax,ic;

  
  char buf[200];
  f1.getline(buf,200);
  
  int i=1;
  
  while(i<=NGRB && f1.getline(buf,200))
    {
      if(sscanf(buf,"%lf %d %lf %lf %lf %lf %lf %lf",&fluence,&nshell,&etot,&r0,&dr0,&gmin,&gmax,&ic)<=0) break;
      i++;
    }
  i--;
  f1.close();

  if(i<NGRB)
    {    
      std::ifstream f2(paramFile.c_str());
      f2.getline(buf,200);
      int nb;
      if((NGRB % i)==0) nb=3;
      else nb=(NGRB % i);
      for(int j = 1; j<=nb;j++)
	{
	  f2.getline(buf,200);
	  sscanf(buf,"%lf %d %lf %lf %lf %lf %lf %lf",&fluence,&nshell,&etot,&r0,&dr0,&gmin,&gmax,&ic);
	}
      f2.close();
    }
  SetParameters(fluence,nshell,etot,r0,dr0,gmin,gmax,ic);
}

void Parameters::SetParameters(double fluence, int nshell, double etot, double r0, double dr0, double gmin, double gmax, double ic)
{
  SetGalDir(-200,-200);
  SetEtot(etot);
  SetFluence(fluence);
  SetInitialSeparation(r0);
  SetInitialThickness(dr0);		
  SetNshell(nshell);
  SetGammaMin(gmin);
  SetGammaMax(gmax);
  SetInverseCompton(ic);
}

void Parameters::ComputeParametersFromFile(std::string paramFile, int NGRB)
{
  using std::pow;
  std::ifstream f1(paramFile.c_str());
  if (!f1.is_open()) 
    {
      std::cout<<"GRBConstants: Error Opening paramFile\n";
      exit(1);
    }
  int    nshell;
  double fluence,r0,dr0;
  double gmin,gmax,etot,ic;
  double tv,ep,Eco;
  
  int GBM;
  
  char buf[200];
  f1.getline(buf,100);
  
   int i=1;
  while(i<=NGRB && f1.getline(buf,100))
    {
      if(sscanf(buf,"%lf %d %lf %lf %lf %lf %d",&fluence,&nshell,&tv,&Eco,&ep,&ic,&GBM)<=0) break;
      i++;
    } 
  i--;
  f1.close();
	
  if(i<NGRB)
    {
      std::ifstream f2(paramFile.c_str());
      f2.getline(buf,100);

      for(int j = 1; j<=(NGRB %i);j++)
	{
	  f2.getline(buf,100);
	  if(sscanf(buf,"%lf %d %lf %lf %lf %lf %d",&fluence,&nshell,&tv,&Eco,&ep,&ic,&GBM)<=0);
	}
      f2.close();
    }
  //////////////////////////////////////////////////
  //  ep/=2.0;
  if (GBM>0) 
    SetGBMOutput(true);
    
  SetGRBNumber(UInt_t(65540+NGRB));  

  if(tv<=0.0)
    {
      if(nshell>1 && nshell<4)
	tv= -2;
      else if(nshell>4)
	tv= -1;
    }
  if(tv==0)   // both
    {
      tv=((rnd->Uniform()<0.3) ? -2 : -1);
      m_Type = 0;
    }
  
  if(tv==-2)
    {
      tv     = pow(10.0,rnd->Gaus(-1.0,0.1)); // short bursts
      m_Type=1;
    }
  else if(tv==-1)
    {
      tv   = pow(10.0,rnd->Gaus(0.0,0.5)); //long burts
      m_Type=2;
    }
  
  //  if(ep==0) SetEpeak(ep);
  //  m_Ep=ep;
  double gmax_gmin=2;  
  double E52   = 1.0;
  double Ep100 = ep/100.0;///pow(gmax_gmin,2.);
  double ae3   = cst::ae * 3.;
  double ab3   = cst::ab * 3.;
  double g100,d7;
  r0    =  2.0 * cst::c * tv;
  double r10 = r0 * 1.0e-10;
  
  double G=0.0;
  
  if(Ep100==0)
    {
      if(Eco==0)
	{
	  g100 = 1.0;
	  d7   = 1.0;
	  G=100.0*g100;
	  dr0=1e7*d7;
	}
      else
	{
	  G = 40.5 * Eco;
	  g100=G/100.0;
	  d7   = 1.0;
	  dr0=1e7*d7;
	}
    }
  else
    {
      if(Eco==0)
	{
	  d7   = 1.0;
	  dr0=1e7*d7;
	  
	  G = 77.8 * pow(E52*ab3/d7,1./4.)*ae3/sqrt(Ep100*tv);
	  g100=G/100.0;
	}
      else
	{
	  G = 40.5 * Eco;
	  d7=13.6*E52*ab3*pow(ae3,4.0)/(pow(Eco,4.0)*pow(Ep100*tv,2.0));
	  dr0=1e7*d7;
	  g100=G/100.0;
	}
    }
  
  Ep100=3.63*sqrt(E52*ab3/d7)*pow(ae3,2.)/(r10*pow(g100,2));
  ep = 100.0*Ep100;
  
  if(DEBUG) 
    std::cout<<"Expected Ep= "<<ep<<" Gamma = "<<G<<" Ecut-off : "<<G/40.5<<" GeV, Rsh = "<<r0*G<<" tv "<<tv<<std::endl;
  ////////////////////////////////////////////////// 
  etot = 1e52*E52;
  gmin = 2.*G/(gmax_gmin + 1.);
  gmax = gmax_gmin*gmin;

  //  SetEpeak(Ep100*100.0);
  //  SetFluence(fluence);
  //  SetGammaMin(gmin);
  //  SetGammaMax(gmax_gmin*gmin);  
  //  SetInitialSeparation(r0);
  //  SetInitialThickness(dr0);
  //  SetNshell(nshell);
  //  SetInverseCompton(ic);
  SetParameters(fluence,nshell,etot,r0,dr0,gmin,gmax,ic);

}

void Parameters::PrintParameters()
{
  std::cout<<"--------------------------------------------------"<<std::endl;
  std::cout<<"-GRB NUMBER :  ---------------> "<<m_GRBnumber<<std::endl;
  std::cout<<" Nshell                      = "<<m_Nshell<<std::endl;
  std::cout<<" Etot                        = "<<m_Etot<<" Erg "<<std::endl;
  std::cout<<" Fluence in the Batse Range  = "<<m_Fluence<<std::endl;
  std::cout<<" Initial Separation          = "<<m_InitialSeparation<<std::endl;
  std::cout<<" Initial Thickness           = "<<m_InitialThickness<<std::endl;
  std::cout<<" Minimum Lorentsz Factor     = "<<m_Gmin<<std::endl;
  std::cout<<" Maximum Lorentsz Factor     = "<<m_Gmax<<std::endl;
  std::cout<<" Inverse Compton Parameter   = "<<m_InverseCompton<<std::endl;
  if(m_GBM) 
    std::cout<<" for this bursts will be generated the GBM output "<<std::endl;
}
