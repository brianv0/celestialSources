#include <vector>
#include <iostream>
#include <fstream>

#include "GRBConstants.h"

Parameters::Parameters()
{
  rnd = new TRandom();
  m_seed = 65539;
  rnd->SetSeed(m_seed);
  m_Nshell=10;
}

void Parameters::SetGRBNumber(UInt_t GRBnumber)
{
  m_seed = GRBnumber;
  rnd->SetSeed(GRBnumber);
}

//////////////////////////////////////////////////
double Parameters::GetBATSEFluence()
{
  if (m_InitialSeparation<pow(10.,9.))
    return pow(10.0,rnd->Gaus(-6.3,0.57)); //erg/cm^2 (Short Bursts)
  return pow(10.0,rnd->Gaus(-5.4,0.62)); //erg/cm^2 (Long Burst)

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
  
  if (nshell<2)
    {
      m_Nshell = int(2+rnd->Integer(100));
      /*
	if (m_InitialSeparation<pow(10,9))
	m_Nshell = (int) pow(10.,rnd->Gaus(1.,0.5)); // 5.0e+8 5.0e8  //short
	else 	
	m_Nshell = (int) pow(10.,rnd->Gaus(2.,0.3)); //5.0e+10 1.0e7 //long
      */
    }
  else
    {
      m_Nshell = nshell;
    }
  
  if(m_Nshell<2 || m_Nshell>500) 
    {
      std::cout<<" Warning setting the shell "<<std::endl;
      SetNshell(0);
    }
}

 void Parameters::SetFluence(double fluence)
{
  if(fluence == 0)  
    m_Fluence = GetBATSEFluence();
  else
    m_Fluence = fluence;
}

void Parameters::SetEtot(double etot)
{
  m_Etot = (etot>0.0) ? etot : pow(10.0,rnd->Gaus(51.,0.5)); //erg
}

void Parameters::SetInitialSeparation(double initialSeparation)
{
  m_InitialSeparation = (initialSeparation>0.0) ? initialSeparation : pow(10.,rnd->Uniform(7.0,10.0));
}

void Parameters::SetInitialThickness(double initialThickness)
{
  //  double logR=log10(m_InitialSeparation);
  m_InitialThickness = (initialThickness>0.0) ? initialThickness : pow(10.,rnd->Uniform(7.0,10.0));
  //pow(10.0,rnd->Uniform(15.0-logR,27.0-2.0*logR));//Gaus(0.0,6.0-0.5*logR));
  //    pow(10.0,15.0-logR+rnd->Uniform(15.0-logR,21.0-1.5*logR));//Gaus(0.0,6.0-0.5*logR));
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
  m_InverseCompton = (ic>=0.0 && ic<=1.0) ? ic : rnd->Uniform(0.,1.);
}

//..................................................

void Parameters::ReadParametersFromFile(std::string paramFile, int NGRB)
{
  
  std::ifstream f1(paramFile.c_str());
  if (!f1.is_open()) 
    {
      std::cout<<"GRBConstants: Error Opening paramFile\n";
      exit(1);
    }
  int    nshell, seed;
  double fluence,etot,r0,dr0,l0,b0;
  double gmin,gmax,ic;

  
  char buf[200];
  f1.getline(buf,200);

  int i=1;

  while(i<=NGRB && f1.getline(buf,200))
    {
      if(sscanf(buf,"%d %lf %lf %lf %d %lf %lf %lf %lf %lf %lf",
		&seed,&l0,&b0,&fluence,&nshell,&etot,&r0,&dr0,&gmin,&gmax,&ic)<=0) break;
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
	  sscanf(buf,"%d %lf %lf %lf %d %lf %lf %lf %lf %lf %lf",
		 &seed,&l0,&b0,&fluence,&nshell,&etot,&r0,&dr0,&gmin,&gmax,&ic);
	}
      f2.close();
    }
  //  SetGRBNumber(seed);
  m_seed = rnd->GetSeed();
  
  SetGalDir(l0,b0);
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
  std::ifstream f1(paramFile.c_str());
  if (!f1.is_open()) 
    {
      std::cout<<"GRBConstants: Error Opening paramFile\n";
      exit(1);
    }
  int    nshell, seed;
  double fluence,etot,r0,dr0,l0,b0;
  double gmin,gmax,ic;
  double gamma,tv,ep;

  char buf[200];
  f1.getline(buf,200);
  
  int i=1;
  
  while(i<=NGRB && f1.getline(buf,200))
    {
      if(sscanf(buf,"%d %lf %lf %lf %d %lf %lf %lf %lf %lf ",
		&seed,&l0,&b0,&fluence,&nshell, &gamma,&tv,&etot,&ep,&ic)<=0) break;
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
	  sscanf(buf,"%d %lf %lf %lf %d %lf %lf %lf %lf %lf ",
		 &seed,&l0,&b0,&fluence,&nshell, &gamma,&tv,&etot,&ep,&ic);
	}
      f2.close();
    }
  //  SetGRBNumber(seed);
  m_seed = rnd->GetSeed();
  
  SetGalDir(l0,b0);
  SetEtot(etot);
  SetFluence(fluence);
  double g100 = gamma/100.;


  r0  = 5.99e10* tv;
  dr0 = 4.00e10 * (etot/1.e52) * pow(tv,-2.) * pow(g100,-4.)* pow(ep,-2.); 

  /*  
  double ec=ep/2.;  
  r0  = 6.e10 * tv;
  dr0 = 8.5e8*pow(etot/1.e52,1./3.)*pow(ec,-1./3.)*pow(ep,-5./3.)*pow(tv,-2./3.);
  double gamma  = 262.0*pow(etot/1.e52,1./6.)*pow(ec/ep,1./12.)*pow(tv,-1/3.);
  */

  gmin = 2./3. * 100.0 * g100;
  SetGammaMin(gmin);
  SetGammaMax(2.*gmin);

  //4.e10 *(etot/1.e52)/pow(ep * tv * pow((gmax-gmin)/100n.,2),2);

  SetInitialSeparation(r0);
  SetInitialThickness(dr0);		
  SetNshell(nshell);
  SetInverseCompton(ic);
  
}



void Parameters::PrintParameters()
{
  std::cout<<"--------------------------------------------------"<<std::endl;
  std::cout<<"-GRB NUMBER :  ---------------> "<<m_seed<<std::endl;
  std::cout<<" Nshell                      = "<<m_Nshell<<std::endl;
  std::cout<<" Etot                        = "<<m_Etot<<" Erg "<<std::endl;
  std::cout<<" Fluence in the Batse Range  = "<<m_Fluence<<std::endl;
  std::cout<<" Initial Separation          = "<<m_InitialSeparation<<std::endl;
  std::cout<<" Initial Thickness           = "<<m_InitialThickness<<std::endl;
  std::cout<<" Minimum Lorentsz Factor     = "<<m_Gmin<<std::endl;
  std::cout<<" Maximum Lorentsz Factor     = "<<m_Gmax<<std::endl;
  std::cout<<" Inverse Compton Parameter   = "<<m_InverseCompton<<std::endl;
}
