#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

#include "GRBConstants.h"

Parameters::Parameters()
{
  rnd = new TRandom();
  SetGRBNumber((long) rnd->GetSeed());
}

void Parameters::SetGRBNumber(UInt_t GRBnumber)
{
  m_GRBnumber = GRBnumber;
  rnd->SetSeed(m_GRBnumber);
  rnd->Uniform();
}

//////////////////////////////////////////////////
double Parameters::GetBATSEFluence()
{
    using std::pow;
  if (m_InitialSeparation<pow(10.,9.))
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
  
  if (nshell<1)
    m_Nshell = int(2+rnd->Integer(100));
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
      rnd->Uniform();  //THB: this needed parentheses
    }
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

  SetGRBNumber(65540+ (long) seed);  

  
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
    using std::pow;
  std::ifstream f1(paramFile.c_str());
  if (!f1.is_open()) 
    {
      std::cout<<"GRBConstants: Error Opening paramFile\n";
      exit(1);
    }
  int    nshell, seed;
  double fluence,r0,dr0,l0,b0;
  double gmin,gmax,ic;
  double gmax_gmin,tv,ep;

  char buf[200];
  f1.getline(buf,100);
  
  int i=1;
  while(i<=NGRB && f1.getline(buf,100))
    {
      if(sscanf(buf,"%d %lf %lf %lf %d %lf %lf %lf %lf %lf ",
		&seed,&l0,&b0,&fluence,&nshell,&tv,&gmax_gmin,&ep,&ic)<=0) break;
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
	  sscanf(buf,"%d %lf %lf %lf %d %lf %lf %lf %lf %lf ",
		 &seed,&l0,&b0,&fluence,&nshell,&tv,&gmax_gmin,&ep,&ic);
	  
	}
      seed=NGRB;
      f2.close();
    }
  if(gmax_gmin<1.5) gmax_gmin=1.5;
  //  if(gmax_gmin>1.5) gmax_gmin=1.5;
  SetGRBNumber(65540+ (long) seed);  
  
  SetGalDir(l0,b0);
  SetEtot(1e52);
  SetFluence(fluence);
  //  double g100 = tau/100.;
  double E52   = 1.0;
  //  ep = pow(10,rnd->Gaus(2.,.5));//00.0)log10(ep),log10(ep/2.0)));
  double Ep100 = ep/100.0/(gmax_gmin*gmax_gmin);
  r0    = cst::c * tv;
  double G100 = 3.0;
  dr0   =  1.647e9*E52 * cst::ab * pow((double)cst::ae,4.0)/(pow(Ep100,2.0)*pow(tv,2.0)*pow(G100,4.0));//*pow(gmax_gmin,2.0);
  //  double gmax_gmin = tau;
  double G = G100*100.0;

  ////////////////////////////////////////////////// 
  /*
  double tc = tv;
  
  r0    = 2. * cst::c * tv;
  dr0   = cst::c * tc;
  
  double G = 484.159 * pow(etot/1.e52,1./4.) * pow(cst::ab,1./4.) * cst::ae * pow(ep,-1./2.) 
    * pow(tv,-1./2.) * pow(tc,-1./4.); 
  double gmax_gmin = 10;
  */
////////////////////////////////////////////////// 
  gmin = 2.*G/(gmax_gmin + 1.);
  //  m_Tau=tau;
  SetGammaMin(gmin);
  SetGammaMax(gmax_gmin*gmin);
  
  SetInitialSeparation(r0);
  SetInitialThickness(dr0);
  SetNshell(nshell);
  SetInverseCompton(ic);
}

/*
  double Parameters::GetNextPeak()
  {
  return rnd->Exp(m_Tau);
  }
*/

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
}
