#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

#include "GRBConstants.h"

#define DEBUG  1

Parameters::Parameters()
  : m_QG(false)
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
  double tmp;
  tmp = rnd->Uniform();
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

double Parameters::GetBATSEDuration()
{
  using std::pow;
  if (m_Type==1)
    return pow(10.0,(double)rnd->Gaus(-0.2,0.55)); // (Short Bursts)
  return pow(10.0,(double)rnd->Gaus(1.46,0.49)); // (Long Burst)
}


//////////////////////////////////////////////////

void Parameters::SetGalDir(double l, double b)
{
  double r1 = rnd->Uniform();
  double r2 = rnd->Uniform();
  
  double ll = (l<=180.0 && l>=-180.0) ? l : 180.-360.*r1;
  double bb = (b<=90.0 && b>=-90.0)   ? b : ((180.0/TMath::Pi())*acos(1.0-2.0*r2)-90.0);
  m_GalDir=std::make_pair(ll,bb);
  m_GalDir=std::make_pair(ll,bb);
}

//////////////////////////////////////////////////

void Parameters::SetFluence(double fluence)
{
  if(fluence == 0)  
    m_Fluence = GetBATSEFluence();
  else
    {
      m_Fluence = fluence;
      double tmp;
      tmp = rnd->Uniform();//THB: this needed parentheses
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

void Parameters::SetParameters(double fluence, double etot, double r0, double dr0, double gmin, double gmax, double ic)
{
  SetGalDir(-200,-200);
  SetEtot(etot);
  SetFluence(fluence);
  SetInitialSeparation(r0);
  SetInitialThickness(dr0);		
  SetGammaMin(gmin);
  SetGammaMax(gmax);
  SetInverseCompton(ic);
}

void Parameters::ComputeParametersFromFile(std::string paramFile, int NGRB)
{
  using std::fabs; using std::pow;
  std::ifstream f1(paramFile.c_str());
  if (!f1.is_open()) 
    {
      std::cout<<"GRBConstants: Error Opening paramFile\n";
      exit(1);
    }
  double fluence,r0,dr0;
  double gmin,gmax,etot,ic;
  double ep,Eco;
  char type;
  
  int GBM,QG;
  
  char buf[200];
  f1.getline(buf,100);
  /// First the parameter file is read
  int i=1;
  while(i<=NGRB && f1.getline(buf,100))
    {
      if(sscanf(buf,"%lf %s %lf %lf %lf %d %d",&fluence,&type,&Eco,&ep,&ic,&GBM,&QG)<=0) break;
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
	  if(sscanf(buf,"%lf %s %lf %lf %lf %d %d",&fluence,&type,&Eco,&ep,&ic,&GBM,&QG)<=0);
	}
      f2.close();
    }
  //////////////////////////////////////////////////
  //  if(QG>0)
  //    std::cout<<"Lorentz Violation effect enabled!"<<std::endl;
  SetQGOutput(QG>0);
  /// Sets the bolean flag for GBM output
  if (GBM>0) 
    SetGBMOutput(true);
  /// Set the type: m_Type=1 for Short bursts, 2 for longs
  /// If type is not decleared, it is extracted randomly
  SetGRBNumber(UInt_t(65540+NGRB));  
  if(type=='S' || type =='s')
    m_Type=1;
  else if(type=='L' || type =='l')
    m_Type=2;
  else 
    m_Type=((rnd->Uniform()<=0.25) ? 1 : 2);
  
  m_Duration = GetBATSEDuration();

  double gmax_gmin=2.0;  
  double E52   = 1.0;
  double Ep100 = ep/100.0;
  double ae3   = cst::ae * 3.;
  double ab3   = cst::ab * 3.;
  double g100=0;
  double r10=0;
  double G=0.0;
  double d7=0.0;
  double tv=0.0;
  
  if (Eco<1.0) Eco = rnd->Uniform(3.0,10.0);
  if(m_Type==1) // SHORT BURSTS
    {
      tv  = m_Duration;
      if(ep==0) ep = pow(10.,log10(235.0)+log10(1.75) * rnd->Gaus());
    }
  else // LONG BURSTS
    {
      tv  = TMath::Max(m_Duration/50.0,pow(10.0,rnd->Gaus(0.0,0.5))); ///FWHM @ 20ke
      if(ep==0) ep = 0.5 * pow(10.,log10(235.0)+log10(1.75) * rnd->Gaus());
    }
  //////////////////////////////////////////////////
  G   =  40.5 * Eco;
  g100=G/100.0;
  
  Ep100 = ep/100.0;
  d7  =  13.6 * E52 * ab3 * pow(ae3,4.0)/(pow(Eco,4.0)*pow(Ep100*tv,2.0)*pow(cst::csi,4.0));
  dr0 =   1e7 * d7;
  
  r0    =  2.0 * cst::c * tv;
  r10   =   r0 * 1.0e-10;
  
  //////////////////////////////////////////////////
  ep = 100.0*Ep100;
  
  if(DEBUG) 
    std::cout<<"Expected Ep= "<<ep<<" Gamma = "<<G<<" Ecut-off : "<<G/40.5<<" GeV, Rsh = "<<r0*G*G<<" tv "<<tv<<std::endl;
  ////////////////////////////////////////////////// 
  etot = 1e52*E52;
  gmin = 2.*G/(gmax_gmin + 1.);
  gmax = gmax_gmin*gmin;

  SetParameters(fluence,etot,r0,dr0,gmin,gmax,ic);

}

void Parameters::PrintParameters()
{
  std::cout<<"--------------------------------------------------"<<std::endl;
  std::cout<<"-GRB NUMBER :  ---------------> "<<m_GRBnumber<<std::endl;
  std::cout<<" Duaration (T90)             = "<<m_Duration<<" s "<<std::endl;
  std::cout<<" Etot                        = "<<m_Etot<<" Erg "<<std::endl;
  std::cout<<" Fluence in the Batse Range  = "<<m_Fluence<<" erg/cm^2/s"<<std::endl;
  std::cout<<" Initial Separation          = "<<m_InitialSeparation<<" cm"<<std::endl;
  std::cout<<" Initial Thickness           = "<<m_InitialThickness<<" cm"<<std::endl;
  std::cout<<" Minimum Lorentsz Factor     = "<<m_Gmin<<std::endl;
  std::cout<<" Maximum Lorentsz Factor     = "<<m_Gmax<<std::endl;
  std::cout<<" Inverse Compton Parameter   = "<<m_InverseCompton<<std::endl;
  if(m_GBM) 
    std::cout<<" for this bursts will be generated the GBM output "<<std::endl;
  if(m_QG)
    std::cout<<" Lorentz Violation set to true for this burst "<<std::endl;
}
