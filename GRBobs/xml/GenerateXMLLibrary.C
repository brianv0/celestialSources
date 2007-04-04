#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>

//////////////////////////////////////////////////
// OPTIONS FOR GENERATING THE XML LIBRARY:
double MinExtractedPhotonEnergy = 10.0; //MeV
long FirstBurstTime  =      1000; //1e4
//double AverageInterval = 86400.0; //[s] 1 grbs / day
//double AverageInterval = 48516.9; //[s] 650 grbs / yr
double AverageInterval = 47351.35135; //[s] 666 grbs / yr

bool  GeneratePF         =   true; // If true: PF is used to normalize Bursts.
                                   // If false Fluence is used to normalize Bursts.
TString GenerateIC        =   "random";//"random"; //"yes", "no"
TString RepointFile       = "none";

double Fssc_Fsyn_const   =   10;
bool  GenerateRedshift   =   true;
bool  GenerateGBM        =   true;//false;
bool  GLASTCoordinate    =   false;
bool  JayDistributions   =   true;
bool  BN                 =   true;
//bool  newXML             =   true;
//////////////////////////////////////////////////
void GenerateGalDir(TRandom3 *rnd, double &l, double &b)
{
  double r1 = rnd->Uniform();
  double r2 = rnd->Uniform();
  l = 180.-360.*r1;
  b = (180.0/TMath::Pi())*acos(1.0-2.0*r2)-90.0;
}



//////////////////////////////////////////////////
class point
{
public:
  point(){;}
  bool is_ARR(){return is_arr;}
  double get_l(){return l;}
  double get_b(){return b;}
  long int get_time(){return time;}

  void set_ARR(bool f){is_arr=f;}
  void set_time(long int  t){time=t;}
  void set_l(double v){l=v;}
  void set_b(double v){b=v;}
  
private:
  bool is_arr;
  long time;
  double l;
  double b;
};



std::vector<point*> getPointingFromFile()
{
  //  std::string filename = "l_b_coord.txt";
  std::ifstream grb_positions(RepointFile,std::ios::in);
  std::vector<point*> positions;
  if(!grb_positions.is_open()) 
    return positions;

  double t, l, b;
  long t0=-1;
      
  while(!grb_positions.eof())
    {
      point *a_position = new point;
      grb_positions >> t;
      grb_positions >> l;
      grb_positions >> b;
      if(long(t)!=t0)
	{
      	  a_position->set_time(long(t));
	  a_position->set_l(l);
	  a_position->set_b(b);
	  a_position->set_ARR(true);
	  positions.push_back(a_position);
	  t0=long(t);
	}
    }
  
  std::cout<<"---------- Pointings: ---------- " <<positions.size()<<std::endl;
  for(unsigned int i=0;i<positions.size();i++)
    std::cout<<i<<" "<<positions[i]->get_time()<<" l= "<<positions[i]->get_l()<<" b= "<<positions[i]->get_b()<<std::endl;
  std::cout<<" -------------------------------------------------- " <<std::endl;
  
  return positions;
}


int GetClosestPointing(std::vector<point*> pointings, long time)
{
  int i=0;
  int N=pointings.size();
  long t = (pointings[i])->get_time();
  while (t<time && i<N-2)
    {
      t = (pointings[++i])->get_time();
    }
  if(time - pointings[i]->get_time() < pointings[i+1]->get_time()-time)
    return i;
  else 
    return i+1;
}

std::vector<point*> GetBurstPositions(int N=1000)
{
  TRandom3 *rnd_2 = new TRandom3();
  rnd_2->SetSeed(65540);
  
  std::vector<point*> pointings        = getPointingFromFile();
  std::vector<point*> bursts_positions;
  
  unsigned S=pointings.size();
  long BurstTime=FirstBurstTime;
  double l, b;
  for (int i=0;i<N;i++)
    {
      point *a_position = new point;
      GenerateGalDir(rnd_2,l,b);
      
      a_position->set_time(BurstTime);
      a_position->set_l(l);
      a_position->set_b(b);
      a_position->set_ARR(false);
      
      BurstTime+=(int) rnd_2->Exp(AverageInterval);
      
      bursts_positions.push_back(a_position);
    }
  for (unsigned i=0;i<S;i++)
    {
      long time = (pointings[i])->get_time();
      int ti = GetClosestPointing(bursts_positions,time);
      bursts_positions[ti]->set_time(pointings[i]->get_time());
      bursts_positions[ti]->set_l(pointings[i]->get_l());
      bursts_positions[ti]->set_b(pointings[i]->get_b());
      bursts_positions[ti]->set_ARR(true);
      //std::cout<<bursts_positions[ti]->get_time()<<", l= "<<bursts_positions[ti]->get_l()<<",b= "<<bursts_positions[ti]->get_b()<<std::endl;
    }
  delete rnd_2;
  for (unsigned i=0;i<S;i++) delete pointings[i];
  std::cout<<" -------------------------------------------------- "<<std::endl;
  std::cout<<"                 Final positions:                   "<<std::endl;
  for(unsigned i=0;i<bursts_positions.size();i++) 
    std::cout<<i<<" "<<bursts_positions[i]->get_time()<<", l= "<<bursts_positions[i]->get_l()<<",b= "<<bursts_positions[i]->get_b()<<std::endl;
  std::cout<<" -------------------------------------------------- "<<std::endl;
  return bursts_positions;
}

double binary1(double z){
  // implemented by D. Noyes, this is a fit to the redshift distribution
  // of simuated binary mergers produced by Chris Fryer.
  
  double val;
  double par1[4] = {1.00690e-01,2.29417e-01,2.04017,8.47315e-01};
  double par2[4] = {-8.76443e-01,4.07510,1.56548e-01,1.21371};
  if(z<=1)
    val = (par1[0]+par1[1]*z+par1[2]*z*z)*exp(-par1[3]*z);
  else
    val = (par2[0]+par2[1]*z+par2[2]*z*z)*exp(-par2[3]*z);
                                                                             
  return val/4.99383068084716797;
                                                                              
}
double GetRedshiftShort(TRandom3 *rnd)
{
  double z=0.0, y=0.0;
  double func=1.;
  do {
    z = rnd->Uniform(0.,10.0);
    y  = rnd->Uniform(0.,1.0);
    func = binary1(z);
  } while (y>func);
 return z;
}


double GRB_SF1a(double z){
  // from Porciani and Madau, ApJ, 526,L522 1999.
  double val;
  val = 0.3*exp(3.4*z)/(exp(3.8*z)+45.)/7.16974794864654541e-01;
  return(val);
}

double GRB_SF2a(double z){
  // from Porciani and Madau, ApJ, 526,L522 1999.
  double val;
  val = 0.15*exp(3.4*z)/(exp(3.4*z)+22.);
  return(val);
}
 
double GRB_SF3a(double z){
  // from Porciani and Madau, ApJ, 526,L522 1999.
  double val;
  val = 0.134*exp(3.05*z)/(exp(2.93*z)+15.);
  val/=2.20879745483398438;
  return(val);
}

double GetRedshiftLong(TRandom3 *rnd)
{
  double z=0.0, y=0.0;
  double func=1.0;
  do {
    z = rnd->Uniform(0.,10.0);
    y  = rnd->Uniform(0.,1.0);
    // you can change the simulated z distrib, by changing GRB_SF1a 
    // to GRB_SF2a or GRB_SF3a.
    func = GRB_SF1a(z);
  } while (y>func);
  return z;
}


double GetPeakFluxLong(TRandom3 *rnd)
{
  // Peak Flux Distributions, from Jay and Jerry IDL code.
  static const int NentriesL = 31;
  //LONG BURSTS
  static const double N[NentriesL] = { 8,      9,     10,     12,     15,     20,     25,     30,     35,   
				       40,     50,     60,     70,     80,     90,    100,    120,    150,   
				       200,    250,    300,    400,    500,    600,    700,    800,    900,   
				       1000,   1100,   1200,   1262 };
  
  static const double P[NentriesL] = { 49.279, 44.093, 42.020, 38.082, 33.853, 28.819, 25.417, 21.994, 19.378,   
				       18.148, 15.910, 13.671, 12.326, 10.387,  9.709,  8.521,  7.118,  5.709,   
				       4.297,  3.366,  2.740,  1.981,  1.503,  1.250,  1.035,  0.858,  0.713,   
				       0.595,  0.500,  0.389,  0.240 };
  
  // Long:
  //  static double MaxP = TMath::MaxElement(NentriesL,P);
  static double MaxN = TMath::MaxElement(NentriesL,N);
  static double MinN = TMath::MinElement(NentriesL,N);
  double Ni = rnd->Uniform(MinN,MaxN);  
  double Ndum=0.0;
  int idum=0;
  
  while(N[idum]<Ni)
    {
      Ndum=N[idum];
      idum++;
    }
  double Fp_hi = P[idum-1];
  double Fp_lo = P[idum];
  double x = rnd->Uniform();
  double   Fp = Fp_lo * pow(1.0 - x * (1.0 - pow(Fp_hi/Fp_lo,-1.0)),-1.0);
  return Fp;
}

double GetPeakFluxShort(TRandom3 *rnd)
{
  // Peak Flux Distributions, from Jay and Jerry IDL code.
  static const int NentriesS = 24;
  //SHORT BURSTS
  static const double M[NentriesS] = { 8,      9,     10,     12,     15,     20,     25,     30,     35,    
				       40,     50,     60,     70,     80,     90,    100,    120,    150,    
				       200,    250,    300,    400,    430,    462 };
  
  static const double Q[NentriesS] = { 27.912, 27.792, 24.635, 21.766, 17.314, 15.150, 13.985, 12.510, 11.102,    
				       9.845,  7.996,  7.318,  6.415,  5.344,  4.795,  4.191,  3.437,  2.998,   
				       2.335,  1.908,  1.572,  1.001,  0.817,  0.371 };
  
  //  static double MaxQ = TMath::MaxElement(NentriesS,Q);
  static double MaxM = TMath::MaxElement(NentriesS,M);
  static double MinM = TMath::MinElement(NentriesS,M);
  double Mi = rnd->Uniform(MinM,MaxM);
  
  double Mdum=0.0;
  int idum=0;
  
  while(M[idum]<Mi)
    {
      Mdum=M[idum];
      idum++;
    }
  double Fp_hi = Q[idum-1];
  double Fp_lo = Q[idum];
  double x = rnd->Uniform();
  double Fp = Fp_lo * pow(1.0 - x * (1.0 - pow(Fp_hi/Fp_lo,-1.0)),-1.0);
  if(Fp<1e-3) std::cout<<" WARNING PF < 1e-3************************************"<<std::endl;
  return Fp;
}

double GetEPeak(TRandom3 *rnd, int Type, double Stretch)
{
  double Ep =  pow(10.,rnd->Gaus(log10(235.0),log10(1.75))); //Short
  if(Type==2 && (GeneratePF)) Ep/=Stretch; //Long
  return Ep;
}
//////////////////////////////////////////////////

void GenerateXMLLibrary(int Nbursts=1000)
{
  std::vector<point*> burst_positions= GetBurstPositions(Nbursts);
  //////////////////////////////////////////////////
  double FL=-1e-5;
  double PF=-1.0;
  double Fluence=FL;
  double PeakFlux=PF;
  double z=0;
  double Duration;
  int BURSTtype = 0; // 1->S, 2->L, else S+L
  // If Glast Coordinae = true:
  double theta = 45.0;
  double phi   = 0.0;
  
  double alpha0 = 10.0;  //  -3<= a < 1.0 
  double beta0  = 10.0;   //   b < a && b < -1 
  double alpha_min = -2.5;  //  -3<= a < 1.0 
  double alpha_max = 0.5;   //   -3<= a < 1.0 
  double beta_max = -1.2;   //   b < a && b < -2 
  double beta_min = -7.2;// -4.0;   //   b < a && b < -2 
  
  //////////////////////////////////////////////////
  int NphLat=0;
  double DelayTime=0;
  double ExtraComponent_Duration=0.0;
  double CO_Energy= 0.0;  
  //////////////////////////////////////////////////
  double Fssc_Fsyn = 0.0;
  double Essc_Esyn = 1e4;
  if(GenerateIC=="yes") 
    Fssc_Fsyn =   Fssc_Fsyn_const;
  //////////////////////////////////////////////////
  if (GenerateRedshift && (CO_Energy!=0)) std::cout<<" WARNING!!! Generate redshift is true and CO_Energy!=0"<<std::endl;
  
  TRandom3 *rnd = new TRandom3();
  rnd->SetSeed(65540);
  
  TRandom3 *rnd_1 = new TRandom3();
  rnd_1->SetSeed(65540);
  
  //  TRandom3 *rnd_2 = new TRandom3();
  //  rnd_2->SetSeed(65540);
  
  //  rnd->SetSeed(19741205);
  //rnd->SetSeed(20030914);
  // Peak Flux Distributions, from Jay and Jerry IDL code.
  
  static const int NentriesL = 31;
  static const int NentriesS = 24;
  
  //LONG BURSTS

  static const double PL[NentriesL] = { 49.279, 44.093, 42.020, 38.082, 33.853, 28.819, 25.417, 21.994, 19.378,   
					18.148, 15.910, 13.671, 12.326, 10.387,  9.709,  8.521,  7.118,  5.709,   
					4.297,  3.366,  2.740,  1.981,  1.503,  1.250,  1.035,  0.858,  0.713,   
					0.595,  0.500,  0.389,  0.240 };
  
  static const double PS[NentriesS] = { 27.912, 27.792, 24.635, 21.766, 17.314, 15.150, 13.985, 12.510, 11.102,    
					9.845,  7.996,  7.318,  6.415,  5.344,  4.795,  4.191,  3.437,  2.998,   
					2.335,  1.908,  1.572,  1.001,  0.817,  0.371 };
  
  double XL[NentriesL];
  double XS[NentriesS];
  
  //////////////////////////////////////////////////
  
  for(int i=0;i<NentriesL;i++)
    {
      XL[i]=PL[NentriesL-1-i];
    }
  
  for(int i=0;i<NentriesS;i++)
    {
      XS[i]=PS[NentriesS-1-i];
    }
  ////////// ALPHA DISTRIBUTION ////////////////////////////////////////
  int pl_histalpha[] = {  1,   2,   3,  11,  11,  27,  33, 100, 130, 143,
			  227, 270, 364, 390, 426, 508, 435, 398, 365, 271,
			  210, 156,  76,  57,  33,  20,  12,  14,  12,   5};
  long sza = sizeof(pl_histalpha)/sizeof(int);
  
  TH1D *AlphaHisto = new TH1D("AlphaHisto","AlphaHisto",sza,-2.5,0.5);  
  for (int i=0; i<sza; i++)
    {
      AlphaHisto->SetBinContent(i+1,pl_histalpha[i]);
    }
  
  ////////// BETA DISTRIBUTION //////////  
  int pl_histbeta[] =  { 18,  18,  18,  18,  18,  18,  18,  18,  18,  18,
			 18,  18,  18,  18,  18,  18,  18,  18,  18,  18,
			 19,  19,  19,  19,  19,  19,  19,  19,  19,  19,
			 19,  19,  20,  20,  35,  37,  27,  40,  75, 122,
			 117, 160, 167, 215, 263, 310, 333, 365, 450, 487,
			 430, 400, 318, 217, 148, 100,  51,  22,   7,   2};
  
  long szb = sizeof(pl_histbeta)/sizeof(int);
  TH1D *BetaHisto = new TH1D("BetaHisto","BetaHisto",szb,-7.2,-1.2);  
  for (int i=0; i<szb; i++)
    {
      BetaHisto->SetBinContent(i+1,pl_histbeta[i]);
    }
  
  //////////////////////////////////////////////////
  

  gDirectory->Delete("PFlong");
  gDirectory->Delete("PFshort");
  
  gDirectory->Delete("FLlong");
  gDirectory->Delete("FLshort");
  
  gDirectory->Delete("T90long");
  gDirectory->Delete("T90short");
  
  gDirectory->Delete("Zlong");
  gDirectory->Delete("Zshort");
  
  gDirectory->Delete("alphalong");
  gDirectory->Delete("alphashort");
  
  gDirectory->Delete("betalong");
  gDirectory->Delete("betashort");
  
  gDirectory->Delete("Eplong");
  gDirectory->Delete("Epshort");

  TH1D* PFlong  = new TH1D("PFlong","",NentriesL-1,XL);
  TH1D* PFshort = new TH1D("PFshort","",NentriesS-1,XS);
  
  TH1D* FLlong  = new TH1D("FLlong","",50,-9,-2);
  TH1D* FLshort = new TH1D("FLshort","",50,-9,-2);
  
  TH1D* T90long  = new TH1D("T90long","",100,-3,3);
  TH1D* T90short = new TH1D("T90short","",100,-3,3);
  
  TH1D* Zlong  = new TH1D("Zlong","",100,0,10);
  TH1D* Zshort = new TH1D("Zshort","",100,0,10);
  
  TH1D* alphalong  = new TH1D("alphalong","",50,alpha_min,alpha_max);
  TH1D* alphashort = new TH1D("alphashort","",50,alpha_min,alpha_max);

  TH1D* betalong  = new TH1D("betalong","",60,beta_min,beta_max);
  TH1D* betashort = new TH1D("betashort","",60,beta_min,beta_max);
  
  TH1D* Eplong  = new TH1D("Eplong","",60,log10(235.0)-5*log10(1.75),log10(235.0)+5*log10(1.75));
  TH1D* Epshort = new TH1D("Epshort","",60,log10(235.0)-5*log10(1.75),log10(235.0)+5*log10(1.75));
  
  TH2D *FluenceVsT90 = new TH2D("FluenceVsT90","",100,-3,3,100,-9,-2);
  TH2D *PeakFluxVsT90 = new TH2D("PeakFluxVsT90","",100,-3,3,100,-2,3);
  
  FluenceVsT90->SetXTitle("Log_{10}(T_{90})");
  FluenceVsT90->SetYTitle("Fluence 50-300 keV (erg/cm^{2}) ");
  
  //////////////////////////////////////////////////
  T90long->GetXaxis()->CenterTitle();
  T90short->GetXaxis()->CenterTitle();

  PFlong->GetXaxis()->CenterTitle();
  PFshort->GetXaxis()->CenterTitle();
  
  betalong->GetXaxis()->CenterTitle();
  betashort->GetXaxis()->CenterTitle();

  alphalong->GetXaxis()->CenterTitle();
  alphashort->GetXaxis()->CenterTitle();
  
  Eplong->GetXaxis()->CenterTitle();
  Epshort->GetXaxis()->CenterTitle();

  Zlong->GetXaxis()->CenterTitle();
  Zshort->GetXaxis()->CenterTitle();

  T90long->GetXaxis()->SetTitle("Log_{10} Duration [sec]");
  T90short->GetXaxis()->SetTitle("Log_{10} Duration [sec]");

  PFlong->GetXaxis()->SetTitle("Peak Flux 50-300 keV [ph/cm^{2}/s]");
  PFshort->GetXaxis()->SetTitle("Peak Flux 50-300 keV [ph/cm^{2}/s]");
  
  betalong->GetXaxis()->SetTitle("High Energy Spectral Index (#beta)");
  betashort->GetXaxis()->SetTitle("High Energy Spectral Index (#beta)");

  alphalong->GetXaxis()->SetTitle("Low Energy Spectral Index (#alpha)");  
  alphashort->GetXaxis()->SetTitle("Low Energy Spectral Index (#alpha)");  
  
  Eplong->GetXaxis()->SetTitle("Peak Of The e^{2} N(e) spectrum (log_{10}) [keV]");
  Epshort->GetXaxis()->SetTitle("Peak Of The e^{2} N(e) spectrum (log_{10}) [keV]");

  Zlong->GetXaxis()->SetTitle("Redshift");
  Zshort->GetXaxis()->SetTitle("Redshift");

  
  if(BN)
    {
      //      PFlong->SetLineStyle(4);
      PFshort->SetLineStyle(3);
      //FLlong->SetLineStyle(4);
      FLshort->SetLineStyle(3);
      //T90long->SetLineStyle(4);
      T90short->SetLineStyle(3);
      //Zlong->SetLineStyle(4);
      Zshort->SetLineStyle(3);
      //alphalong->SetLineStyle(4);
      alphashort->SetLineStyle(3);
      //betalong->SetLineStyle(4);
      betashort->SetLineStyle(3);
      //Eplong->SetLineStyle(4);
      Epshort->SetLineStyle(3);
    }
  else
    {

  PFlong->SetLineColor(4);
  PFshort->SetLineColor(2);
  FLlong->SetLineColor(4);
  FLshort->SetLineColor(2);
  T90long->SetLineColor(4);
  T90short->SetLineColor(2);
  Zlong->SetLineColor(4);
  Zshort->SetLineColor(2);
  alphalong->SetLineColor(4);
  alphashort->SetLineColor(2);
  betalong->SetLineColor(4);
  betashort->SetLineColor(2);
  Eplong->SetLineColor(4);
  Epshort->SetLineColor(2);
    }
  Zlong->SetXTitle("Redshift z");
  Zshort->SetXTitle("Redshift z");

  std::ofstream osXML("GRBobs_user_library.xml",std::ios::out);
  std::ofstream osTest("../src/test/GRBParam.txt",std::ios::out);
  osXML.precision(4);
  
  osXML<<"<!-- $Header -->"<<std::endl;
  osXML<<"<!-- ************************************************************************** -->"<<std::endl;
  osXML<<"<source_library title=\"GRBobs_user_library\">"<<std::endl;
  
  long BurstTime=FirstBurstTime;
  //  FILE osTest("../src/test/GRBParam.txt");
  
  
  if(GeneratePF)
    {
      //    fprint(osTest," T0  T90   Z    PF  alpha   beta   Eco ");
      osTest<<setw(11)<<"T0"<<setw(9)<<"T90"<<setw(9)<<"PF"<<setw(9)<<"z"<<setw(9)<<"alpha"<<setw(9)<<"beta"<<setw(10)<<"Ep(keV)";
      osTest<<setw(9)<<" Essc/Esyn"<<setw(9)<<" Fssc/Fsyn "<<setw(9)<<"Eco(GeV)"<<std::endl;
    }
  else
    {
      osTest<<setw(11)<<"T0"<<setw(9)<<"T90"<<setw(9)<<"FL"<<setw(9)<<"z"<<setw(9)<<"alpha"<<setw(9)<<"beta"<<setw(10)<<"Ep(keV)";
      osTest<<setw(9)<<" Essc/Esyn"<<setw(9)<<" Fssc/Fsyn "<<setw(9)<<"Eco(GeV)"<<std::endl;
      //osTest<<" T0  T90    FL  alpha   beta   Eco "<<std::endl;
    }
  
  int type;
  int Nlong=0;
  int Nshort=0;
  
  // From jay et al:
  double nfracL = 0.7;
  double log_mean_long  = log10(15.0);
  double log_mean_short = log10(0.4);
  double log_sigma_long = 0.6;
  double log_sigma_short = 0.4;
  //////////////////////////////////////////////////
  const double a = 0.277;
  const double b = -0.722;
  const double c = 1.47;


  for(int i = 0; i<Nbursts ; i++)
    {
      double gal_l = 0.0; 
      double gal_b = 0.0; 
      //GenerateGalDir(rnd_2,gal_l,gal_b);
      gal_l = burst_positions[i]->get_l();
      gal_b = burst_positions[i]->get_b();
      BurstTime = burst_positions[i]->get_time();
      
      if(BURSTtype==1) 
	{
	  type=1;
	  rnd->Uniform();
	}
      else if(BURSTtype==2) 
	{
	  type=2;
	  rnd->Uniform();
	}
      else type = ((rnd->Uniform() > nfracL) ? 1 : 2);

      //////////////////////////////////////////////////
      theta = 0.0;//10.0*(i%8);
      phi   = 0.0;//10.0*(i%37);
      double alpha = alpha0;
      double beta  = beta0;
      double Ep    = 1.0;
      
      // Handle the Autonomous repoint case.
      double my_beta_min = beta_min;
      if (burst_positions[i]->is_ARR())
	{
	  my_beta_min= -2.0;
	  
	  if(GenerateIC=="random")
	    {
	      double randomNumber = rnd_1->Uniform();
	      if(randomNumber      > 0.6)
		Fssc_Fsyn =   10.0;
	    }
	  else 
	    Fssc_Fsyn = 0.0;
	}
      else if(GenerateIC=="random")
	{
	  double randomNumber = rnd_1->Uniform();
	  if(randomNumber > 0.95)
	    Fssc_Fsyn =    1.0;
	  else
	    Fssc_Fsyn =    0.0;
	}

      
      if(JayDistributions)
	{
	  while (alpha < alpha_min || alpha > alpha_max) alpha = AlphaHisto->GetRandom();
	  while(beta >= alpha || beta >= beta_max || beta < my_beta_min) beta = BetaHisto->GetRandom();
	}
      else
	{
	  while (alpha < alpha_min || alpha > alpha_max) alpha = rnd->Gaus(-1.0,0.4);
	  while(beta >= alpha || beta >= beta_max || beta < my_beta_min) beta = rnd->Gaus(-2.25,0.4);
	  
	}
      
      double PFRND_S       = GetPeakFluxShort(rnd); //ph/cm^2/s (Short Bursts)
      double FluenceRND_S  = pow(10.0,(double)rnd->Gaus(-6.8,0.43));
      double Duration_S    = pow(10.0,(double)rnd->Gaus(log_mean_short,log_sigma_short)); //(Short Bursts)
      double zRND_S        = GetRedshiftShort(rnd);
      
      double Stretch;
      double PFRND_L       = GetPeakFluxLong(rnd); //ph/cm^2/s (Short Bursts)
      double FluenceRND_L  = pow(10.0,(double)rnd->Gaus(-5.5,0.61)); //erg/cm^2 (Long Burst)
      double Duration_L    = 1e6;
      while(Duration_L>640.0)
	{
	  Duration_L = pow(10.0,(double)rnd->Gaus(log_mean_long,log_sigma_long)); //(Long Bursts)
	}
      double zRND_L        = GetRedshiftLong(rnd);
      //////////////////////////////////////////////////
      double co_energy = CO_Energy;
      if(CO_Energy < 0.0) co_energy = floor(rnd->Uniform(1.,10.));
      
      //////////////////////////////////////////////////
      if (type==1) //SHORT
	{
	  //// Peak flux normalization:
	  //	  long seed = rnd->GetSeed();
	  if(GeneratePF)
	    {
	      if(PF<=0) PeakFlux = PFRND_S; //ph/cm^2/s (Short Bursts)
	      else PeakFlux=PF; //ph/cm^2/s (Short Bursts)
	      PFshort->Fill(PeakFlux);
	    }
	  else
	    {
	      if(FL<=0)
		Fluence = FluenceRND_S;
	      else
		Fluence = FL; //erg/cm^2 (Short Bursts)
	      FLshort->Fill(log10(Fluence));
	      if(Fluence>1e-3) std::cout<<" WARNING F > 1e-3************************************"<<std::endl;
	    }
	  
	  if (GenerateRedshift)
	    z = zRND_S;
	  else
	    z = 0.0;
	  
	  Duration = Duration_S;
	  Ep = GetEPeak(rnd, type, 1.0);
	  
	  T90short->Fill(log10(Duration));
	  Zshort->Fill(z);
	  alphashort->Fill(alpha);
	  betashort->Fill(beta); 
	  Epshort->Fill(log10(Ep));
	  Nshort++;
	}
      else //LONG
	{
	  if(GeneratePF)
	    {	
	      if(PF<=0)
		PeakFlux = PFRND_L;
	      else
		PeakFlux =  PF; //ph/cm^2/s (Long Bursts)
	      
	      double lpf = log10(PeakFlux);
	      Stretch =  TMath::Max(1.0, a * lpf*lpf + b * lpf + c);
	      PFlong->Fill(PeakFlux);
	    }
	  else
	    {
	      if(FL<=0)
		Fluence = FluenceRND_L;
	      else
		Fluence = FL;
	      Stretch =  1.0;
	      FLlong->Fill(log10(Fluence));
	      if(Fluence>1e-3) std::cout<<" WARNING F > 1e-3************************************"<<std::endl;
	    }
	  
	  
	  Duration = Duration_L * Stretch;
	  T90long->Fill(log10(Duration));
	  
	  
	  if (GenerateRedshift)
	    z = zRND_L;
	  else
	    z = 0.0;
	  
	  
	  Ep=GetEPeak(rnd, type, Stretch);
	  
	  alphalong->Fill(alpha);
	  betalong->Fill(beta); 
	  Eplong->Fill(log10(Ep));
	  Zlong->Fill(z);
	  Nlong++;
	}
      
      if(GeneratePF) 
	PeakFluxVsT90->Fill(log10(Duration),log10(PeakFlux));
      else 
	FluenceVsT90->Fill(log10(Duration),log10(Fluence));
      
      osXML<<""<<std::endl;
      if(i<10) osXML<<"<source name=\" GRB_0000"<<i<<" \">"<<std::endl;
      else if(i<100) osXML<<"<source name=\" GRB_000"<<i<<" \">"<<std::endl;
      else if(i<1000) osXML<<"<source name=\" GRB_00"<<i<<" \">"<<std::endl;
      else if(i<10000) osXML<<"<source name=\" GRB_0"<<i<<" \">"<<std::endl;
      else osXML<<"<source name=\" GRB_"<<i<<" \">"<<std::endl;
      osXML<<"<spectrum escale=\"MeV\">"<<std::endl;
      
      
      if(GeneratePF)
	osXML<<" <SpectrumClass name=\"GRBobsmanager\" params=\"tstart="<<BurstTime<<", duration="<<Duration<<", peakFlux="<<PeakFlux;
      else 
	osXML<<" <SpectrumClass name=\"GRBobsmanager\" params=\"tstart="<<BurstTime<<", duration="<<Duration<<", fluence="<<Fluence;
      
      //      if(!GLASTCoordinate)  
	osXML<<", l="<<gal_l<<", b="<<gal_b;
      
      osXML<<", redshift="<<z<<", alpha="<<alpha<<", beta="<<beta<<", Ep="<<Ep<<", emin="<<MinExtractedPhotonEnergy;
      
      if(Fssc_Fsyn>0.0)
	osXML<<", essc_esyn="<<Essc_Esyn<<", Fssc_Fsyn="<<Fssc_Fsyn;
      osXML<<", GBM="<<(int)GenerateGBM;
      if(NphLat>0)
	osXML<<", EC_NLAT="<<NphLat<<", EC_delay="<<DelayTime<<", EC_duration="<<ExtraComponent_Duration;
      if (co_energy>0)
	osXML<<", EC_CutOff="<<co_energy;
      osXML<<" \"/>"<<std::endl;
      
      if(GLASTCoordinate)  
	osXML<<"<direction theta=\""<<theta<<"\" phi=\""<<phi<<"\" />"<<std::endl; //u galactic_dir l=\""<<l<<"\" b=\""<<b<<"\" />"<<std::endl;
      else
	osXML<<"<use_spectrum frame=\"galaxy\" />"<<std::endl; //u galactic_dir l=\""<<l<<"\" b=\""<<b<<"\" />"<<std::endl;
      
      osXML<<"</spectrum> </source>"<<std::endl;
      //////////////////////////////////////////////////
      
      if(GeneratePF) 
	osTest<<setw(11)<<BurstTime<<" "<<setw(8)<<setprecision(4)<<Duration<<" "<<setw(8)<<PeakFlux<<" ";
      else
	osTest<<setw(11)<<BurstTime<<" "<<setw(8)<<setprecision(4)<<Duration<<" "<<setw(8)<<Fluence<<" ";
      
      osTest<<setw(8)<<z<<" "<<setw(8)<<alpha<<" "<<setw(8)<<beta<<" "<<setw(8)<<Ep<<" "<<setw(8)<<Essc_Esyn<<setw(8)<<Fssc_Fsyn<<setw(8)<<co_energy<<std::endl;
      //////////////////////////////////////////////////
      //      BurstTime+=(int) rnd->Exp(AverageInterval);
      //      time = GetClosesTime(pointings,BurstTime);
    
      
    }
  //////////////////////////////////////////////////
  osXML<<" "<<std::endl;
  osXML<<"<source name=\" GRBobs-10\" >"<<std::endl;
  for(int i = 0; i<10 ; i++)
    {
      if(i<10) osXML<<"     <nestedSource sourceRef=\" GRB_0000"<<i<<" \"/>"<<std::endl;
    }
  
  osXML<<"</source>"<<std::endl;
  //////////////////////////////////////////////////
  
  //////////////////////////////////////////////////
  osXML<<""<<std::endl;
  osXML<<"<source name=\" GRBobs-100\" >"<<std::endl;
  for(int i = 0; i<100 ; i++)
    {
      if(i<10) osXML<<"     <nestedSource sourceRef=\" GRB_0000"<<i<<" \"/>"<<std::endl;
      else if(i<100) osXML<<"     <nestedSource sourceRef=\" GRB_000"<<i<<" \"/>"<<std::endl;
      else osXML<<"     <nestedSource sourceRef=\" GRB_"<<i<<" \"/>"<<std::endl;
    }
  osXML<<"</source>"<<std::endl;
  //////////////////////////////////////////////////
  
  //////////////////////////////////////////////////
  osXML<<""<<std::endl;
  osXML<<"<source name=\" GRBobs-1000\" >"<<std::endl;
  for(int i = 0; i<1000 ; i++)
    {
      if(i<10) osXML<<"     <nestedSource sourceRef=\" GRB_0000"<<i<<" \"/>"<<std::endl;
      else if(i<100) osXML<<"     <nestedSource sourceRef=\" GRB_000"<<i<<" \"/>"<<std::endl;
      else if(i<1000) osXML<<"     <nestedSource sourceRef=\" GRB_00"<<i<<" \"/>"<<std::endl;
    }
  osXML<<"</source>"<<std::endl;
  //////////////////////////////////////////////////
  
  osXML<<""<<std::endl;
  osXML<<"<source name=\" AllGRBobs\" >"<<std::endl;
  for(int i = 0; i<Nbursts ; i++)
    {
      if(i<10) osXML<<"     <nestedSource sourceRef=\" GRB_0000"<<i<<" \"/>"<<std::endl;
      else if(i<100) osXML<<"     <nestedSource sourceRef=\" GRB_000"<<i<<" \"/>"<<std::endl;
      else if(i<1000) osXML<<"     <nestedSource sourceRef=\" GRB_00"<<i<<" \"/>"<<std::endl;
      else if(i<10000) osXML<<"     <nestedSource sourceRef=\" GRB_0"<<i<<" \"/>"<<std::endl;
      else osXML<<"     <nestedSource sourceRef=\" GRB_"<<i<<" \"/>"<<std::endl;
    }
  osXML<<"</source>"<<std::endl;
  osXML<<"  </source_library>"<<std::endl;
  TCanvas *Correlation = new TCanvas("Correlation","Correlation",600,400);
  if(GeneratePF) PeakFluxVsT90->Draw();
  else  FluenceVsT90->Draw();
    
  /*
    TCanvas *Redshifts = new TCanvas("Redshifts","Redshifts",600,400);
    Redshifts->SetFillColor(10);
  if(Nlong>Nshort){
    Zlong->Draw();
    Zshort->Draw("same");
  }else{
    Zshort->Draw();
    Zlong->Draw("same");
  }
  */
  TCanvas *Distributions = new TCanvas("Distributions","Distributions",900,600);
  Distributions->SetFillColor(10);
  gStyle->SetOptStat(0);
  Distributions->Divide(3,2);
  if(Nlong>Nshort)
    {
      Distributions->cd(1);
      T90long->Draw();
      T90short->Draw("same");
      
      Distributions->cd(2);
      if(GeneratePF)
	{
	  gPad->SetLogx();
	  gPad->SetLogy();
	  PFlong->Draw();
	  PFshort->Draw("same");
	}
      else
	{
	  FLlong->Draw();
	  FLshort->Draw("same");
	}
      Distributions->cd(3);
      Zlong->Draw();
      Zshort->Draw("same");

      Distributions->cd(4);
      alphalong->Draw();
      alphashort->Draw("same");
      
      Distributions->cd(5);
      betalong->Draw();
      betashort->Draw("same");
      
      Distributions->cd(6);
      Eplong->Draw();
      Epshort->Draw("same");
    }
  else
    {
      Distributions->cd(1);
      
      T90short->Draw();
      T90long->Draw("same");
      
      Distributions->cd(2);
      gPad->SetLogx();
      gPad->SetLogy();
      PFshort->Draw();
      PFlong->Draw("same");

      Distributions->cd(3);
      Zshort->Draw();
      Zlong->Draw("same");

      Distributions->cd(4);
      alphashort->Draw();
      alphalong->Draw("same");
      
      Distributions->cd(5);
      betashort->Draw();
      betalong->Draw("same");
      
      Distributions->cd(6);
      Epshort->Draw();
      Eplong->Draw("same");
    }
  std::cout<<"--------------------------------------------------"<<std::endl;
  std::cout<<"Generate "<<Nlong<<" Long Burst"<<std::endl;
  std::cout<<"Generate "<<Nshort<<" Short Burst"<<std::endl;
  std::cout<<"--------------------------------------------------"<<std::endl;
  
  Distributions->cd();
  Distributions->Print("GeneratedDistributions.eps");
  Correlation->Print("GeneratedCorrelations.eps");
  //  if (GenerateRedshift)  Redshifts->Print("GeneratedRedshift.eps");

  //  for (int i=0;i<burst_positions.size();i++) delete burst_positions[i];
}



