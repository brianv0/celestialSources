#include <iomanip>

double GetRedshiftLong(TRandom *rnd)
{
  double z=0.0;
  while(z<=0) z = rnd->Gaus(2.0,1.0);
  return z;
}

double GetRedshiftShort(TRandom *rnd)
{
  double z=0.0;
  while(z<=0) z = rnd->Gaus(.1,.1);
  return z;
}

double GetPeakFluxLong(TRandom *rnd)
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
  static double MaxP = TMath::MaxElement(NentriesL,P);
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

double GetPeakFluxShort(TRandom *rnd)
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

  static double MaxQ = TMath::MaxElement(NentriesS,Q);
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

//////////////////////////////////////////////////
void GenerateXMLLibrary(int Nbursts=100)
{
  TRandom *rnd = new TRandom();
  rnd->SetSeed(65540);
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
  

  double MinExtractedPhotonEnergy = 30.0; //MeV
  long FirstBurstTime  =      1000; //1e4
  double AverageInterval = 86400.0; //s
  bool  GeneratePF  =   false;//true;// If true: PF is used to normalize Bursts.
  double FL=-1e-5;
  double PF=-1.0;
  double Fluence,PeakFlux;
  double z;
  double Duration;
  int BURSTtype = 0; // 1->S, 2->L, else S+L
  bool GLASTCoordinate            = false;
  double theta = 45.0;
  double phi   = 0.0;
  double alpha0 = 10.0;  //  -3<= a < 1.0 
  double beta0  = 10.0;   //   b < a && b < -1 
  double alpha_min = -3.0;  //  -3<= a < 1.0 
  double alpha_max = 1.0;   //   -3<= a < 1.0 
  double beta_max = -1.5;   //   b < a && b < -2 
  double beta_min = -7.2;   //   b < a && b < -2 
  int  GBM                        = 1;
  
  int NphLat=0;
  double DelayTime=0;
  double ExtraComponent_Duration=0.0;
  double CO_Energy= 0.0;  
  

  gDirectory->Delete("PFlong");
  gDirectory->Delete("PFshort");

  gDirectory->Delete("FLlong");
  gDirectory->Delete("FLshort");

  gDirectory->Delete("T90long");
  gDirectory->Delete("T90short");

  TH1D* PFlong  = new TH1D("PFlong","Peak Flux ",NentriesL-1,XL);
  TH1D* PFshort = new TH1D("PFshort","Peak Flux",NentriesS-1,XS);
  
  TH1D* FLlong  = new TH1D("FLlong","Fluence ",50,-9,-2);
  TH1D* FLshort = new TH1D("FLshort","Fluence",50,-9,-2);

  TH1D* T90long  = new TH1D("T90long","Duration",100,-3,3);
  TH1D* T90short = new TH1D("T90short","Duration",100,-3,3);

  TH1D* alphalong  = new TH1D("alphalong","Low energy spectral index",100,alpha_min,alpha_max);
  TH1D* alphashort = new TH1D("alphashort","Low energy spectral index",100,alpha_min,alpha_max);
  TH1D* betalong  = new TH1D("betalong","High energy spectral index",100,-4.0,beta_max);
  TH1D* betashort = new TH1D("betashort","High energy spectral index",100,-4.0,beta_max);
  
  PFlong->SetLineColor(4);
  PFshort->SetLineColor(2);
  FLlong->SetLineColor(4);
  FLshort->SetLineColor(2);
  T90long->SetLineColor(4);
  T90short->SetLineColor(2);
  alphalong->SetLineColor(4);
  alphashort->SetLineColor(2);
  betalong->SetLineColor(4);
  betashort->SetLineColor(2);

  
  std::ofstream os("GRBobs_user_library.xml",std::ios::out);
  std::ofstream osTest("../src/test/GRBParam.txt",std::ios::out);
  os.precision(4);

  os<<"<!-- $Header -->"<<std::endl;
  os<<"<!-- ************************************************************************** -->"<<std::endl;
  os<<"<source_library title=\"GRBobs_user_library\">"<<std::endl;
  char SourceName[100]="GRB";
  
  long BurstTime=FirstBurstTime;

  //  FILE osTest("../src/test/GRBParam.txt");
  

  if(GeneratePF)
    //    fprint(osTest," T0  T90   Z    PF  alpha   beta   Eco ");
    osTest<<setw(11)<<"T0"<<setw(9)<<"T90"<<setw(9)<<"PF"<<setw(9)<<"z"<<setw(9)<<"alpha"<<setw(9)<<"beta"<<setw(10)<<"Eco"<<std::endl;
  //osTest<<" T0  T90    PF  alpha   beta   Eco "<<std::endl;
  else
    //    fprint(osTest," T0  T90   Z    PF  alpha   beta   Eco ");
    osTest<<setw(11)<<"T0"<<setw(9)<<"T90"<<setw(9)<<"FL"<<setw(9)<<"z"<<setw(9)<<"alpha"<<setw(9)<<"beta"<<setw(10)<<"Eco"<<std::endl;
  //osTest<<" T0  T90    FL  alpha   beta   Eco "<<std::endl;
  
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
      while (alpha < alpha_min || alpha > alpha_max) alpha = rnd->Gaus(-1.0,0.4);
      while(beta >= alpha || beta >= beta_max || beta < beta_min) beta = rnd->Gaus(-2.25,0.4);
      //////////////////////////////////////////////////
      if (type==1) //SHORT
	{
	  //// Peak flux normalization:
	  long seed = rnd->GetSeed();
	  if(GeneratePF)
	    {
	      if(PF<=0)
		PeakFlux = GetPeakFluxShort(rnd); //ph/cm^2/s (Short Bursts)
	      else
		{
		  PeakFlux=PF; //ph/cm^2/s (Short Bursts)
		  GetPeakFluxShort(rnd);
		}
	      PFshort->Fill(PeakFlux);	  
	      //	  Duration = pow(10.0,(double)rnd->Gaus(-0.2,0.55)); //(Short Bursts)
	    }
	  else
	    {
	      if(FL<=0)
		{
		  Fluence = pow(10.0,(double)rnd->Gaus(-6.3,0.57));
		}
	      else
		{
		  Fluence = FL; //erg/cm^2 (Short Bursts)
		  rnd->Gaus();
		}
	      FLshort->Fill(log10(Fluence));
	      if(Fluence>1e-3) std::cout<<" WARNING F > 1e-3************************************"<<std::endl;
	    }

	  rnd->SetSeed(seed);
	  Duration = pow(10.0,(double)rnd->Gaus(log_mean_short,log_sigma_short)); //(Short Bursts)
	  z = GetRedshiftShort(rnd);

	  T90short->Fill(log10(Duration));
	  alphashort->Fill(alpha);
	  betashort->Fill(beta); 
	  Nshort++;
	}
      else //LONG
	{
	  double Stretch;
	  long seed = rnd->GetSeed();
	  if(GeneratePF)
	    {	
	      if(PF<=0)
		{
		  PeakFlux = GetPeakFluxLong(rnd);
		}
	      else
		{
		  PeakFlux =  PF; //ph/cm^2/s (Long Bursts)
		}
	      double lpf = log10(PeakFlux);
	      Stretch =  TMath::Max(1.0, a * lpf*lpf + b * lpf + c);
	      PFlong->Fill(PeakFlux);
	    }
	  else
	    {
	      if(FL<=0)
		{
		  Fluence = pow(10.0,(double)rnd->Gaus(-5.4,0.62)); //erg/cm^2 (Long Burst)
		}
	      else
		{
		  Fluence = FL;
		}
	      Stretch =  1.0;
	      FLlong->Fill(log10(Fluence));
	      if(Fluence>1e-3) std::cout<<" WARNING F > 1e-3************************************"<<std::endl;
	    }
	  
	  Duration=1e6;
	  rnd->SetSeed(seed);
	  while(Duration>640.0)
	    {
	      Duration = pow(10.0,(double)rnd->Gaus(log_mean_long,log_sigma_long)); //(Long Bursts)
	      Duration*=Stretch;
	    }
	  T90long->Fill(log10(Duration));
	  z = GetRedshiftLong(rnd);
	  alphalong->Fill(alpha);
	  betalong->Fill(beta); 
	  
	  Nlong++;

	}
      
      os<<""<<std::endl;
      if(i<10) os<<"<source name=\" GRB_0000"<<i<<" \">"<<std::endl;
      else if(i<100) os<<"<source name=\" GRB_000"<<i<<" \">"<<std::endl;
      else if(i<1000) os<<"<source name=\" GRB_00"<<i<<" \">"<<std::endl;
      else if(i<10000) os<<"<source name=\" GRB_0"<<i<<" \">"<<std::endl;
      else os<<"<source name=\" GRB_"<<i<<" \">"<<std::endl;
      os<<"<spectrum escale=\"MeV\">"<<std::endl;

      double co_energy = CO_Energy;
      if(co_energy < 0.0) co_energy = floor(rnd->Uniform(1.,10.));
      
      if(GeneratePF)
	os<<" <SpectrumClass name=\"GRBobsmanager\" params=\""<<BurstTime<<" , "<<Duration<<" , "<<PeakFlux<<" , "<<
	  z<<" , "<<alpha<<" , "<<beta<<" , "<<MinExtractedPhotonEnergy<<" , "<<GBM;
      else 
	os<<" <SpectrumClass name=\"GRBobsmanager\" params=\""<<BurstTime<<" , "<<Duration<<" , "<<Fluence<<" , "<<
	  z<<" , "<<alpha<<" , "<<beta<<" , "<<MinExtractedPhotonEnergy<<" , "<<GBM;
      
      os<<" , "<<NphLat<<" , "<<DelayTime<<" , "<<ExtraComponent_Duration<<" , "<<co_energy<<" \"/>"<<std::endl;
      



      if(GLASTCoordinate)  
	os<<"<direction theta=\""<<theta<<"\" phi=\""<<phi<<"\" />"<<std::endl;//u galactic_dir l=\""<<l<<"\" b=\""<<b<<"\" />"<<std::endl;
      else
	os<<"<use_spectrum frame=\"galaxy\" />"<<std::endl;//u galactic_dir l=\""<<l<<"\" b=\""<<b<<"\" />"<<std::endl;
      
      
      os<<"</spectrum> </source>"<<std::endl;
      //////////////////////////////////////////////////
      
      
      if(GeneratePF) 
	osTest<<setw(11)<<BurstTime<<" "<<setw(8)<<setprecision(4)<<Duration<<" "<<setw(8)<<PeakFlux<<" "<<setw(8)<<z<<" "<<setw(8)<<alpha<<" "<<setw(8)<<beta<<" "<<setw(8)<<co_energy<<std::endl;
      else
	osTest<<setw(11)<<BurstTime<<" "<<setw(8)<<setprecision(5)<<Duration<<" "<<setw(8)<<Fluence<<" "<<setw(8)<<z<<" "<<setw(8)<<alpha<<" "<<setw(8)<<beta<<" "<<setw(8)<<co_energy<<std::endl;
      //////////////////////////////////////////////////
      
      
      /*
	if(GeneratePF) 
	fprintf("%d %.3f  %.3f %.3f %.3f %.3f %d",BurstTime,Duration,PeakFlux,z,alpha,beta,co_energy);
	else
	fprintf("%d %.3f  %.3f %.3f %.3f %.3f %d",BurstTime,Duration,Fluence,z,alpha,beta,co_energy);
      */
      //////////////////////////////////////////////////
      BurstTime+=(int) rnd->Exp(AverageInterval);
    }
  //////////////////////////////////////////////////
  os<<""<<std::endl;
  os<<"<source name=\" GRBobs-10\" >"<<std::endl;
  for(int i = 0; i<10 ; i++)
    {
      if(i<10) os<<"     <nestedSource sourceRef=\" GRB_0000"<<i<<" \"/>"<<std::endl;
    }
  os<<"</source>"<<std::endl;
  //////////////////////////////////////////////////
  
  //////////////////////////////////////////////////
  os<<""<<std::endl;
  os<<"<source name=\" GRBobs-100\" >"<<std::endl;
  for(int i = 0; i<100 ; i++)
    {
      if(i<10) os<<"     <nestedSource sourceRef=\" GRB_0000"<<i<<" \"/>"<<std::endl;
      else if(i<100) os<<"     <nestedSource sourceRef=\" GRB_000"<<i<<" \"/>"<<std::endl;
      else os<<"     <nestedSource sourceRef=\" GRB_"<<i<<" \"/>"<<std::endl;
    }
  os<<"</source>"<<std::endl;
  //////////////////////////////////////////////////
  
  //////////////////////////////////////////////////
  os<<""<<std::endl;
  os<<"<source name=\" GRBobs-1000\" >"<<std::endl;
  for(int i = 0; i<1000 ; i++)
    {
      if(i<10) os<<"     <nestedSource sourceRef=\" GRB_0000"<<i<<" \"/>"<<std::endl;
      else if(i<100) os<<"     <nestedSource sourceRef=\" GRB_000"<<i<<" \"/>"<<std::endl;
      else if(i<1000) os<<"     <nestedSource sourceRef=\" GRB_00"<<i<<" \"/>"<<std::endl;
    }
  os<<"</source>"<<std::endl;
  //////////////////////////////////////////////////
  
  os<<""<<std::endl;
  os<<"<source name=\" AllGRBobs\" >"<<std::endl;
  for(int i = 0; i<Nbursts ; i++)
    {
      if(i<10) os<<"     <nestedSource sourceRef=\" GRB_0000"<<i<<" \"/>"<<std::endl;
      else if(i<100) os<<"     <nestedSource sourceRef=\" GRB_000"<<i<<" \"/>"<<std::endl;
      else if(i<1000) os<<"     <nestedSource sourceRef=\" GRB_00"<<i<<" \"/>"<<std::endl;
      else if(i<10000) os<<"     <nestedSource sourceRef=\" GRB_0"<<i<<" \"/>"<<std::endl;
      else os<<"     <nestedSource sourceRef=\" GRB_"<<i<<" \"/>"<<std::endl;
    }
  os<<"</source>"<<std::endl;
  os<<"  </source_library>"<<std::endl;

  
  TCanvas *Distributions = new TCanvas("Distributions","Distributions",500,800);
  Distributions->SetFillColor(10);
  gStyle->SetOptStat(0);
  Distributions->Divide(2,2);
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
      alphalong->Draw();
      alphashort->Draw("same");
      
      Distributions->cd(4);
      betalong->Draw();
      betashort->Draw("same");
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
      alphashort->Draw();
      alphalong->Draw("same");
      
      Distributions->cd(4);
      betashort->Draw();
      betalong->Draw("same");
    }
  std::cout<<"--------------------------------------------------"<<std::endl;
  std::cout<<"Generate "<<Nlong<<" Long Burst"<<std::endl;
  std::cout<<"Generate "<<Nshort<<" Short Burst"<<std::endl;
  std::cout<<"--------------------------------------------------"<<std::endl;

  Distributions->cd();
  Distributions->Print("GeneratedDistributions.eps");


}



