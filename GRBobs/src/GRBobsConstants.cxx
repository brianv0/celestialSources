#include <vector>
#include <iostream>
#include <fstream>

#include "../GRBobs/GRBobsConstants.h"
using namespace ObsCst;

GRBobsParameters::GRBobsParameters()
{
  rnd = new TRandom();
  SetGRBNumber((long) rnd->GetSeed());
  m_enph=emin;
}

//////////////////////////////////////////////////
double GRBobsParameters::GetBATSEFluence()
{
  return pow(10.0,rnd->Gaus(-5.4,0.62)); //erg/cm^2 (Long Burst)
  //return pow(10.0,rnd->Gaus(-6.3,0.57); //erg/cm^2 (Short Bursts)
}
//////////////////////////////////////////////////

void GRBobsParameters::SetGRBNumber(long GRBnumber)
{
  m_GRBnumber = GRBnumber;
  rnd->SetSeed(m_GRBnumber);
  rnd->Uniform();
}

void GRBobsParameters::SetFluence(double fluence)
{
  m_fluence = fluence;
  if(m_fluence<=0)  m_fluence = GetBATSEFluence();
}


void GRBobsParameters::SetNumberOfPulses(int NumberOfPulses)
{
  if(NumberOfPulses>0) 
    {
      m_numberOfPulses = NumberOfPulses;
    }
  else 
    {
      double type = rnd->Uniform();
      m_numberOfPulses = (type>0.3) ? (int) pow(10.0,rnd->Gaus(1.47,0.53)) : 1;
    }
}

void GRBobsParameters::SetMinPhotonEnergy(double enph)
{
  m_enph = TMath::Min(emax, TMath::Max(enph,emin));
}


//////////////////////////////////////////////////
void GRBobsParameters::SetGalDir(double l, double b)
{
  
  double ll = (l<=180.0 && l>=-180.0) ? l : rnd->Uniform(-180.0,180.0);
  double bb = (b<=90.0 && b>=-90.0)   ? b : rnd->Uniform(-90.0,90.0);
  m_GalDir=std::make_pair(ll,bb);
}

void GRBobsParameters::GenerateParameters()
{
  m_Peakedness      = pow(10.0,rnd->Gaus(0.16,0.3));
  m_FWHM            = pow(10.0,rnd->Gaus(-0.1,0.5)); //FWHM @ 20keV 
  m_pulseSeparation = pow(10.0,rnd->Gaus(0.13,0.4));
  m_decayTime       = 0.75*pow(0.69,1./m_Peakedness)*m_FWHM;
  m_riseTime        = 0.33*m_decayTime;
  m_pulseHeight     = rnd->Uniform();
  m_Epeak           = pow(10.,rnd->Gaus(2.5,0.1));

}

void GRBobsParameters::PrintParameters()
{
  std::cout<<" N = "<<m_numberOfPulses<<" f = "<<m_fluence<<" dt = "<<m_decayTime<<" rt = "<<m_riseTime
	   <<" pe = "<<m_Peakedness<<" ph = "<<m_pulseHeight<<" ps = "<<m_pulseSeparation
	   <<" ep = "<<m_Epeak<<" a = "<<m_LowEnergy<<" b = "<<m_HighEnergy<<std::endl;
}




//..................................................//

void GRBobsParameters::ReadParametersFromFile(std::string paramFile, int NGRB)
{
  
  std::ifstream f1(paramFile.c_str());
  if (!f1.is_open()) 
    {
      std::cout<<"GRBobsConstants: Error Opening paramFile\n";
      exit(1);
    }
  
  int seed,NumberOfPulses;
  double l0,b0,fluence,alpha,beta;
  
  char buf[100];
  f1.getline(buf,100);
  
  int i=1;
  
  while(i<=NGRB && f1.getline(buf,100))
    {
      if(sscanf(buf,"%d %lf %lf %lf %d %lf %lf",
		&seed,&l0,&b0,&fluence,&NumberOfPulses,&alpha,&beta)<=0) break;
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
	  sscanf(buf,"%d %lf %lf %lf %d %lf %lf",
		 &seed,&l0,&b0,&fluence,&NumberOfPulses,&alpha,&beta);
	}
      seed=NGRB;
      f2.close();
    }
  SetGRBNumber(65540+ (long) seed);  
  SetGalDir(l0,b0);
  SetNumberOfPulses(NumberOfPulses);
  SetFluence(fluence);
  SetAlphaBeta(alpha,beta);
}
