//////////////////////////////////////////////////
void GenerateXMLLibrary(int Nbursts=100)
{
  // Peak Flux Distributions, from Jay and Jarry IDL code.
  static const int NentriesL = 31;
  static const int NentriesS = 24;
  //LONG BURSTS
  double N[NentriesL] = { 8,      9,     10,     12,     15,     20,     25,     30,     35,   
			  40,     50,     60,     70,     80,     90,    100,    120,    150,   
			  200,    250,    300,    400,    500,    600,    700,    800,    900,   
			  1000,   1100,   1200,   1262 };
  
  double P[NentriesL] = { 49.279, 44.093, 42.020, 38.082, 33.853, 28.819, 25.417, 21.994, 19.378,   
			  18.148, 15.910, 13.671, 12.326, 10.387,  9.709,  8.521,  7.118,  5.709,   
			  4.297,  3.366,  2.740,  1.981,  1.503,  1.250,  1.035,  0.858,  0.713,   
			  0.595,  0.500,  0.389,  0.240 };
  
  //SHORT BURSTS
  double M[NentriesS] = { 8,      9,     10,     12,     15,     20,     25,     30,     35,    
			  40,     50,     60,     70,     80,     90,    100,    120,    150,    
			  200,    250,    300,    400,    430,    462 };
  
  double Q[NentriesS] = { 27.912, 27.792, 24.635, 21.766, 17.314, 15.150, 13.985, 12.510, 11.102,    
			  9.845,  7.996,  7.318,  6.415,  5.344,  4.795,  4.191,  3.437,  2.998,   
			  2.335,  1.908,  1.572,  1.001,  0.817,  0.371 };
  
  double P1[NentriesL];
  double Q1[NentriesS];
  
  for(int i=0;i<NentriesL;i++)
    {
      P1[i]=P[NentriesL-1-i];
    }
  for(int i=0;i<NentriesS;i++)
    {
      Q1[i]=Q[NentriesS-1-i];
    }
  
  gDirectory->Delete("PFlong");
  gDirectory->Delete("PFshort");
  
  TH1D* PFlong  = new TH1D("PFlong","PFlong",NentriesL-1,P1);
  TH1D* PFshort = new TH1D("PFshort","PFshort",NentriesS-1,Q1);
  
  for(int i=1; i<NentriesL; i++)
    {
      PFlong->SetBinContent(i,N[NentriesL-i]-N[NentriesL-1-i]);
    }
  PFlong->SetBinContent(NentriesL,0.0);
  for(int i=1;i<NentriesS;i++)
    {
      PFshort->SetBinContent(i+1,M[NentriesS-i]-M[NentriesS-i-1]);
    }
  PFshort->SetBinContent(NentriesS,0.0);
  //////////////////////////////////////////////////
  
  
  
  double MinExtractedPhotonEnergy = 30.0; //MeV
  double FirstBurstTime  =      1000; //1e4
  double AverageInterval = 2000;//86400.0; //s
  bool  GenerateFluence  =   true;//false;
  bool  GeneratePF  =   true;//false;
  double FL=1e-5;
  double PF=1.0;
  double Fluence,PeakFlux;
  double Duration;
  int BURSTtype = 0; // 1->S, 2->L, else S+L
  bool GLASTCoordinate            = false;
  double theta = 45.0;
  double phi   = 0.0;
  
  int  GBM                        = 0;
  
  std::ofstream os("GRBobs_user_library.xml",std::ios::out);
  std::ofstream osTest("../src/test/GRBParam.txt",std::ios::out);
  
  os<<"<!-- $Header -->"<<std::endl;
  os<<"<!-- ************************************************************************** -->"<<std::endl;
  os<<"<source_library title=\"GRBobs_user_library\">"<<std::endl;
  char SourceName[100]="GRB";
  
  double BurstTime=FirstBurstTime;
  
  double alpha;
  double beta;
  if(GeneratePF)
    osTest<<" T0  PF T90  alpha beta"<<std::endl;
  else 
    osTest<<" T0  FL T90  alpha beta"<<std::endl;
  
  int type;
  int Nlong=0;
  int Nshort=0;
  
  // From jay et al:
  double nfracL = 0.7;
  double log_mean_long  = log10(15.0);
  double log_mean_short = log10(0.4);
  double log_sigma_long = 0.6;
  double log_sigma_short = 0.4;
  TRandom *rnd = new TRandom();
  for(int i = 0; i<Nbursts ; i++)
    {
      if(BURSTtype==1) type=1;
      else if(BURSTtype==2) type=2;
      else type = ((rnd->Uniform() > nfracL) ? 1 : 2);
      
      if (type==1)
	{
	  Duration = pow(10.0,(double)rnd->Gaus(log_mean_short,log_sigma_short)); //(Short Bursts)
	  //	  Duration = pow(10.0,(double)rnd->Gaus(-0.2,0.55)); //(Short Bursts)
	  Fluence  = ((GenerateFluence) ? pow(10.0,(double)rnd->Gaus(-6.3,0.57)) : FL); //erg/cm^2 (Short Bursts)
	  PeakFlux = ((GeneratePF)      ? PFlong->GetRandom(): PF); //ph/cm^2/s (Short Bursts)
	  Nshort++;
	}
      else
	{
	  Duration = pow(10.0,(double)rnd->Gaus(log_mean_long,log_sigma_long)); //(Long Bursts)
	  //	  Duration = pow(10.0,(double)rnd->Gaus(1.46,0.49)); //(Long Burst)
	  Fluence = ((GenerateFluence) ? pow(10.0,(double)rnd->Gaus(-5.4,0.62)): FL); //erg/cm^2 (Long Burst)
	  PeakFlux = ((GeneratePF)      ? PFshort->GetRandom(): PF); //ph/cm^2/s (Long Bursts)
	  Nlong++;
	}
      
      theta = 0.0;//10.0*(i%8);
      phi   = 0.0;//10.0*(i%37);
      alpha= rnd->Uniform(0.0, 1.5);
      beta = -1.0*rnd->Uniform(2.0, 2.5);
      os<<""<<std::endl;
      if(i<10) os<<"<source name=\" GRB_0000"<<i<<" \">"<<std::endl;
      else if(i<100) os<<"<source name=\" GRB_000"<<i<<" \">"<<std::endl;
      else if(i<1000) os<<"<source name=\" GRB_00"<<i<<" \">"<<std::endl;
      else if(i<10000) os<<"<source name=\" GRB_0"<<i<<" \">"<<std::endl;
      else os<<"<source name=\" GRB_"<<i<<" \">"<<std::endl;
      os<<"<spectrum escale=\"MeV\">"<<std::endl;
      if(GeneratePF)
	os<<" <SpectrumClass name=\"GRBobsmanager\" params=\""<<BurstTime<<" , "<<Duration<<" , "<<PeakFlux<<" , "<<alpha<<" , "<<beta<<" , "<<MinExtractedPhotonEnergy<<" , "<<GBM<<" \"/>"<<std::endl;
      else 
	os<<" <SpectrumClass name=\"GRBobsmanager\" params=\""<<BurstTime<<" , "<<Duration<<" , "<<Fluence<<" , "<<alpha<<" , "<<beta<<" , "<<MinExtractedPhotonEnergy<<" , "<<GBM<<" \"/>"<<std::endl;

      if(GLASTCoordinate)  
	os<<"<direction theta=\""<<theta<<"\" phi=\""<<phi<<"\" />"<<std::endl;//u galactic_dir l=\""<<l<<"\" b=\""<<b<<"\" />"<<std::endl;
      else
	os<<"<use_spectrum frame=\"galaxy\" />"<<std::endl;//u galactic_dir l=\""<<l<<"\" b=\""<<b<<"\" />"<<std::endl;
      os<<"</spectrum> </source>"<<std::endl;
      //////////////////////////////////////////////////
      if(GeneratePF) 
	osTest<<BurstTime<<" "<<Duration<<" "<<PeakFlux<<" "<<alpha<<" "<<beta<<std::endl;
      else
	osTest<<BurstTime<<" "<<Duration<<" "<<Fluence<<" "<<alpha<<" "<<beta<<std::endl;
      //////////////////////////////////////////////////
      BurstTime+=rnd->Exp(AverageInterval);
    }
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
  
  std::cout<<"--------------------------------------------------"<<std::endl;
  std::cout<<"Generate "<<Nlong<<" Long Burst"<<std::endl;
  std::cout<<"Generate "<<Nshort<<" Short Burst"<<std::endl;
  std::cout<<"--------------------------------------------------"<<std::endl;

}



