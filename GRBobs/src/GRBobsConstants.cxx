#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

#include "../GRBobs/GRBobsConstants.h"

using namespace ObsCst;
using std::pow;

GRBobsParameters::GRBobsParameters()
{
  rnd = new TRandom();
  SetGRBNumber((long) rnd->GetSeed());
  m_Type=0; //1->Short, 2->Long, 0->Both
  m_enph=emin;
}
//////////////////////////////////////////////////
// Relevant BATSE distributions
double GRBobsParameters::GetBATSEDuration()
{
  if (m_Type==1)
    return pow(10.0,(double)rnd->Gaus(-0.2,0.55)); //erg/cm^2 (Short Bursts)
  return pow(10.0,(double)rnd->Gaus(1.46,0.49)); //erg/cm^2 (Long Burst)
}

double GRBobsParameters::GetBATSEFluence()
{
  if (m_Type==1)
    return pow(10.0,(double)rnd->Gaus(-6.3,0.57)); //erg/cm^2 (Short Bursts)
  return pow(10.0,(double)rnd->Gaus(-5.4,0.62)); //erg/cm^2 (Long Burst)
}
//////////////////////////////////////////////////

void GRBobsParameters::SetGRBNumber(long GRBnumber)
{
  m_GRBnumber = GRBnumber;
  rnd->SetSeed(m_GRBnumber);
  double tmp;
  tmp = rnd->Uniform();
}

void GRBobsParameters::SetFluence(double fluence)
{
  m_fluence = fluence;
  if(m_fluence<=0)  m_fluence = GetBATSEFluence();
}


void GRBobsParameters::SetDuration(double duration)
{
  m_duration = duration;
  if(m_duration<2.0) m_Type=1;
  else m_Type=2;
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
  m_RD=0.0;
  m_Peakedness      = 1.5;//pow(10.0,rnd->Gaus(0.16,0.3));
  m_FWHM            = (m_Type==1) ? rnd->Uniform()*m_duration : pow(10.0,rnd->Gaus(-0.1,0.5)); //FWHM @ 20keV 
  m_pulseSeparation = pow(10.0,rnd->Gaus(0.13,0.4)); // double distribution, to be done...
  while(m_RD<=0) m_RD = rnd->Gaus(0.4,0.1);
  m_decayTime       = 1.00/(1.0+m_RD)*pow(0.69,-1./m_Peakedness)*m_FWHM;
  m_riseTime        = m_RD/(1.0+m_RD)*pow(0.69,-1./m_Peakedness)*m_FWHM;
  m_pulseHeight     = rnd->Uniform();
  m_Epeak           = pow(10.,rnd->Gaus(log10(235.0),log10(1.75))); //Short
  if(m_Type==2) m_Epeak/=2.0; //Long
}

void GRBobsParameters::PrintParameters()
{
  std::cout<<" Parameters: Duration = "<<m_duration<<" f = "<<m_fluence<<" dt = "<<m_decayTime<<" rt = "<<m_riseTime
	   <<" pe = "<<m_Peakedness<<" ph = "<<m_pulseHeight<<" tau = "<<m_pulseSeparation
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
  double tstart;
  double duration;
  double fluence,alpha,beta;
  
  char buf[100];
  f1.getline(buf,100);
  
  int i=1;
  
  while(i<=NGRB && f1.getline(buf,100))
    {
      if(sscanf(buf,"%lf %lf %lf %lf %lf",&tstart,&duration,&fluence,&alpha,&beta)<=0) break;
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
	  sscanf(buf,"%lf %lf %lf %lf %lf",&tstart,&duration,&fluence,&alpha,&beta);
	}
      f2.close();
    }
  SetGRBNumber(65540+ (long) floor(tstart));
  SetDuration(duration); // This determines the type.
  SetFluence(fluence);
  SetAlphaBeta(alpha,beta);
  
  SetMinPhotonEnergy(3e4); //keV (this is a defaul value)
  SetGalDir(-200,-200);
  SetGRBNumber(65540+ (long) floor(tstart));
}
