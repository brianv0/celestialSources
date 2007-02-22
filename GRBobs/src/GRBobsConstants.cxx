#include "../GRBobs/GRBobsConstants.h"

using namespace ObsCst;
using std::pow;

GRBobsParameters::GRBobsParameters()
{
  rnd = new TRandom();
  rndGalacticDir = new TRandom();
  SetGRBNumber((long) rnd->GetSeed());
  m_Type=0; //1->Short, 2->Long, 0->Both
  m_enph=enph;
  m_NormType='P';
  m_Stretch=1.0;
}

//////////////////////////////////////////////////
void GRBobsParameters::SetDuration(double duration)
{
  m_duration = duration;
  if(m_duration<=2) 
    m_Type = 1;
  else 
    {
      m_Type = 2;
      //      if(NormType=='P')	
      //	m_duration*=m_Stretch;
    }

}

void GRBobsParameters::SetFluence(double fluence)
{
  if (fluence>1.0e-3) 
    {
      SetPeakFlux(fluence);
      m_NormType='P';
    }
  else 
    {
      m_NormType='F';
      m_fluence = fluence;
    }
}

void GRBobsParameters::SetPeakFlux(double peakflux)
{
  m_peakFlux = peakflux;
  double lpf = log10(m_peakFlux);
  double a = 0.277;
  double b = -0.722;
  double c = 1.47;
  m_Stretch =  TMath::Max(1.0, a * lpf*lpf + b * lpf + c);
}

void GRBobsParameters::SetGRBNumber(long GRBnumber)
{
  m_GRBnumber = GRBnumber;
  rnd->SetSeed(m_GRBnumber);
  rndGalacticDir->SetSeed(m_GRBnumber);
  double tmp;
  tmp = rnd->Uniform();
  tmp = rndGalacticDir->Uniform();
  
}

void GRBobsParameters::SetMinPhotonEnergy(double enph)
{
  m_enph = enph;
}


//////////////////////////////////////////////////
void GRBobsParameters::SetGalDir(double l, double b)
{
  double r1 = rndGalacticDir->Uniform();
  double r2 = rndGalacticDir->Uniform();
  double ll = (l<=180.0 && l>=-180.0) ? l : 180.-360.*r1;
  double bb = (b<=90.0 && b>=-90.0)   ? b : ((180.0/TMath::Pi())*acos(1.0-2.0*r2)-90.0);
  m_GalDir=std::make_pair(ll,bb);
}

double GRBobsParameters::GetBATSEFWHM()
{
  double minwid = 0.01;
  double maxwid = 10.;
  double logfac0 = pow(maxwid / minwid,1./14.);
  double *LogWidth = new double[15];
  double Ngtwid[15]= {430,427,424,421,417,410,388,334,248,178,119, 81, 46, 15,  2};
  for(int i=0;i<15;i++)
    {
      LogWidth[i]= minwid * pow(logfac0,i);      
    }
  double maxNgtwid = TMath::MaxElement(15,Ngtwid);
  double minNgtwid = TMath::MinElement(15,Ngtwid);
  double pickN= rnd->Uniform(minNgtwid,maxNgtwid);
  int idx=0;
  while(Ngtwid[idx] > pickN) idx++;
  
  double wid_lo = LogWidth[idx-1];
  double wid_hi = LogWidth[idx];
  double alpha;
  if (wid_lo <= 0.10)
    alpha = -0.1;
  else if (wid_lo <= 0.25)
    alpha = -1.5;
  else if (wid_lo <= 0.40)
    alpha = 0.25;
  else alpha =  0.6;
  delete[] LogWidth;
  return wid_lo * pow(1. - rnd->Uniform() *(1. - pow(wid_hi/wid_lo,-alpha)),-1./alpha);
}

void GRBobsParameters::GenerateParameters()
{
  m_RD                = RD_Ratio;

  while(m_RD<=0) 
    m_RD = rnd->Gaus(0.4,0.1);
  
  m_Peakedness        = Peakedness;
  while (Peakedness==0) 
    m_Peakedness = pow(10.0,rnd->Gaus(0.16,0.3));
  
  if(rnd->Uniform()<episode_pulses)
    m_pulseSeparation = pow(10.0,rnd->Gaus(log10(pulse_mean_interval),logsigma));
  else
    m_pulseSeparation = pow(10.0,rnd->Gaus(log10(episode_mean_interval),logsigma));
  
  m_FWHM =  GetBATSEFWHM();
  
  if(m_Type==1) m_FWHM/=10.0;
  
  m_decayTime         = 1.00/(1.0+m_RD) / pow(log10(2.0),1./m_Peakedness) * m_FWHM;
  m_riseTime          = m_RD/(1.0+m_RD) / pow(log10(2.0),1./m_Peakedness) * m_FWHM;
  m_pulseHeight       = rnd->Uniform();
  if(m_Epeak == 0 )
    {  
      m_Epeak             =  pow(10.,rnd->Gaus(log10(235.0),log10(1.75))); //Short
      if(m_Type==2 && m_NormType=='P') m_Epeak/=m_Stretch; //Long
    }
}

void GRBobsParameters::PrintParameters()
{
  std::cout<<" Parameters: Duration = "<<m_duration;
  if(m_NormType=='P')  
    std::cout<<" PF = "<<m_peakFlux;
  else 
    std::cout<<" FL = "<<m_fluence;

  std::cout<<" dt = "<<m_decayTime<<" rt = "<<m_riseTime
	   <<" pe = "<<m_Peakedness<<" ph = "<<m_pulseHeight<<" tau = "<<m_pulseSeparation
	   <<" ep = "<<m_Epeak<<" a = "<<m_LowEnergy<<" b = "<<m_HighEnergy<<std::endl;
  std::cout<<"Redshift = "<<m_z<<std::endl;
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
  double fluence;  
  double z;
  double alpha;
  double beta;
  double Ep;
  double fssc_fsyn,essc_esyn;
  double Eco;

  char buf[100];
  f1.getline(buf,100);
  
  int i=1;
  
  while(i<=NGRB && f1.getline(buf,100))
    {
      if(sscanf(buf,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&tstart,&duration,&fluence,&z,&alpha,&beta,&Ep,&essc_esyn,&fssc_fsyn, &Eco)<=0) break;
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
	  sscanf(buf,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&tstart,&duration,&fluence,&z,&alpha,&beta, &Ep,&essc_esyn,&fssc_fsyn, &Eco);
	}
      f2.close();
    }
  
  SetGRBNumber(65540+ (long) floor(tstart));
  
  SetFluence(fluence);
  //  SetPeakFlux(fluence);
  SetDuration(duration);
  SetAlphaBeta(alpha,beta);
  SetEpeak(Ep);
  SetFssc_Fsyn(fssc_fsyn);
  SetEssc_Esyn(essc_esyn);
  SetMinPhotonEnergy(3e4); //keV (this is a defaul value)
  SetGalDir(-200,-200);
  SetCutOffEnergy(Eco);
  SetRedshift(z);
  SetGRBNumber(65540+ (long) floor(tstart));
  
}
