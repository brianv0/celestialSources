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
void Parameters::SetGalDir(double l, double b)
{
  double ll,bb;
  
  ll = (l<=180.0 && l>=-180.0) ? l : rnd->Uniform(180.0,-180.0);
  bb = (b<=180.0 && b>=-180.0) ? b : rnd->Uniform(90.0,-90.0);
  
  m_GalDir=std::make_pair(ll,bb);
}

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

int Parameters::ReadParametersFromFile(std::string paramFile, int NGRB)
{
  
  ifstream f1(paramFile.c_str());
  if (!f1.is_open()) 
    {
      std::cout<<"Error Opening paramFile\n";
      exit(1);
    }
  int    nshell, seed;
  double fluence,etot,r0,dr0,l0,b0;
  double gmin,gmax,ic;


  char buf[100];
  f1.getline(buf,100);
  /*
    printf(buf);
    std::cout<<" "<<std::endl;
    std::cout<<"--------------------------------------------------"<<std::endl;
  */  
  int i=1;

  while(i<NGRB && f1.getline(buf,100))
    {
      i++;
    }
  i--;
  //  SetGRBNumber((long)GetGRBNumber()+1);
  
  f1>>seed>>l0>>b0>>fluence>>nshell>>etot>>r0>>dr0>>gmin>>gmax>>ic;
  //cout<<nshell<<endl;
  if(nshell>1000)
    {
      f1.close();
      f1.open(paramFile.c_str());

      f1.getline(buf,100);
      /*
	printf(buf);
	cout<<" "<<endl;
      */
      //      cout<<"NGRB= "<<NGRB<<", i = "<<i<<" NGRB % i = "<<(NGRB % i)<<endl;
      for(int j = 0; j< (NGRB %i);j++)
	{
	  f1.getline(buf,100);
	}
      
      f1>>seed>>l0>>b0>>fluence>>nshell>>etot>>r0>>dr0>>gmin>>gmax>>ic;    
    }
  
  /*
    cout<<seed<<" "<<l0<<" "<<b0
    <<" "<<fluence<<" "<<nshell<<" "<<
    etot<<" "<<r0<<" "<<dr0<<" "<<gmin<<" "<<gmax<<" "<<ic<<endl;
  */
  
  SetGalDir(l0,b0);
  SetNshell(nshell);
  SetEtot(etot);
  SetFluence(fluence);
  SetInitialSeparation(r0);
  SetInitialThickness(dr0);		
  SetGammaMin(gmin);
  SetGammaMax(gmax);
  SetInverseCompton(ic);

  if(seed>0) SetGRBNumber(65540+ (long) seed);
  
}

/*
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
  double fluence,etot,r0,dr0;
  double gmin,gmax,ic;

  f1.getline(buf,50);
  sscanf(buf,"%d",&nshell);

  f1.getline(buf,50);
  sscanf(buf,"%lf",&etot);

  f1.getline(buf,50);
  sscanf(buf,"%lf",&fluence);
  
  f1.getline(buf,50);
  sscanf(buf,"%lf",&r0);

  f1.getline(buf,50);
  sscanf(buf,"%lf",&dr0);

  
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
  SetInitialSeparation(r0);
  SetInitialThickness(dr0);		
  SetGammaMin(gmin);
  SetGammaMax(gmax);
  SetInverseCompton(ic);
  if(seed>0) SetGRBNumber(GetGRBNumber()+ (long) seed);

}
*/

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
