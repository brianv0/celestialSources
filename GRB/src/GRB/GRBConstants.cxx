#include "GRBConstants.h"

Parameters::Parameters()
{
  rnd = new TRandom();
  SetGRBNumber((long) rnd->GetSeed());
}

double Parameters::GetBATSEFluence()
{
  return pow(10.0,rnd->Gaus(-6.00,0.5)); //erg/cm^2
}

void Parameters::SetNshell(int nshell)
{
  m_nshell = (nshell>0) ? nshell : 10;
}
void Parameters::SetFluence(double fluence)
{
  m_fluence = (fluence>0.0) ? fluence : GetBATSEFluence();
}

void Parameters::SetEtot(double etot)
{
  m_etot = (etot>0.0) ? etot : pow(10.0,rnd->Gaus(50.,53.0)); //erg
}


void Parameters::SetInitialSeparation(double initialSeparation)
{
  m_initialSeparation = (initialSeparation>0.0) ? initialSeparation : 1.0e10;
}
void Parameters::SetInitialThickness(double initialThickness)
{
  m_initialThickness = (initialThickness>0.0) ? initialThickness : 1.0e10;
}

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

  int    nshell;
  double fluence,etot,initialSeparation,initialThickness;  
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
  //cout<<nshell<<" "<<etot<<" "<<fluence<<" "<<initialSeparation<<" "<<initialThickness<<endl;
  SetNshell(nshell);
  SetEtot(etot);
  SetFluence(fluence);
  SetInitialSeparation(initialSeparation);
  SetInitialThickness(initialThickness);		
}

void Parameters::PrintParameters()
{
  std::cout<<"--------------------------------------------------"<<std::endl;
  std::cout<<" Nshell                      = "<<m_nshell<<std::endl;
  std::cout<<" Etot                        = "<<m_etot<<" Erg "<<std::endl;
  std::cout<<" Fluence in the Batse Range  = "<<m_fluence<<std::endl;
  std::cout<<" Initial Separation          = "<<m_initialSeparation<<std::endl;
  std::cout<<" Initial Thickness          = "<<m_initialThickness<<std::endl;
}
