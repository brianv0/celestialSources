#include "GRBConstants.h"

Parameters::Parameters()
{
  rnd = new TRandom();
  SetGRBNumber((long) rnd->GetSeed());
}
//////////////////////////////////////////////////
double Parameters::GetBATSEFluence()
{
  return pow(10.0,rnd->Gaus(-6.00,0.5)); //erg/cm^2
}

double Parameters::GetLESI()
{
  //alpha 
  return rnd->Gaus(-1.0,0.5); 
}

double Parameters::GetHESI()
{
  // beta 
  return rnd->Gaus(-2.5,0.5); 
}
//////////////////////////////////////////////////
void Parameters::SetNshell(int nshell)
{
  m_Nshell = (nshell>0) ? nshell : 10;
}

void Parameters::SetFluence(double fluence)
{
  m_Fluence = (fluence>0.0) ? fluence : GetBATSEFluence();
}

void Parameters::SetEtot(double etot)
{
  m_Etot = (etot>0.0) ? etot : pow(10.0,rnd->Gaus(50.,53.0)); //erg
}

void Parameters::SetInitialSeparation(double initialSeparation)
{
  m_InitialSeparation = (initialSeparation>0.0) ? initialSeparation : 1.0e10;
}

void Parameters::SetInitialThickness(double initialThickness)
{
  m_InitialThickness = (initialThickness>0.0) ? initialThickness : 1.0e10;
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
void Parameters::ReadParametersFromFile(std::string paramFile)
{
  char buf[50];
  ifstream f1(paramFile.c_str());
  if (! f1.is_open()) 
    {
      std::cout<<"Error Opening paramFile.c_str()\n";
      exit(1);
    }
  SetGRBNumber((long)GetGRBNumber()+1);

  int    nshell, seed;
  double fluence,etot,initialSeparation,initialThickness;
  double gmin,gmax,ic;

  f1.getline(buf,50);
  sscanf(buf,"%d",&nshell);

  f1.getline(buf,50);
  sscanf(buf,"%lf",&etot);

  f1.getline(buf,50);
  sscanf(buf,"%lf",&fluence);
  
  f1.getline(buf,50);
  sscanf(buf,"%lf",&initialSeparation);

  f1.getline(buf,50);
  sscanf(buf,"%lf",&initialThickness);

  
  f1.getline(buf,50);
  sscanf(buf,"%lf",&gmin);
  
  f1.getline(buf,50);
  sscanf(buf,"%lf",&gmax);

  f1.getline(buf,50);
  sscanf(buf,"%lf",&ic);

  f1.getline(buf,5);
  sscanf(buf,"%d",&seed);

  SetNshell(nshell);
  SetEtot(etot);
  SetFluence(fluence);
  SetInitialSeparation(initialSeparation);
  SetInitialThickness(initialThickness);		
  SetGammaMin(gmin);
  SetGammaMax(gmax);
  SetInverseCompton(ic);
  if(seed>0) SetGRBNumber(GetGRBNumber()+ (long) seed);

}

void Parameters::PrintParameters()
{
  std::cout<<"--------------------------------------------------"<<std::endl;
  std::cout<<" Nshell                      = "<<m_Nshell<<std::endl;
  std::cout<<" Etot                        = "<<m_Etot<<" Erg "<<std::endl;
  std::cout<<" Fluence in the Batse Range  = "<<m_Fluence<<std::endl;
  std::cout<<" Initial Separation          = "<<m_InitialSeparation<<std::endl;
  std::cout<<" Initial Thickness           = "<<m_InitialThickness<<std::endl;
  std::cout<<" Minimum Lorentsz Factor     = "<<m_Gmin<<std::endl;
  std::cout<<" Maximum Lorentsz Factor     = "<<m_Gmax<<std::endl;
  std::cout<<" Inverse Compton Parameter   = "<<m_InverseCompton<<std::endl;
}
