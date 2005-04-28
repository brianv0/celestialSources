//////////////////////////////////////////////////
void GenerateXMLLibrary(int N=100)
{
  TRandom *rnd = new TRandom();
  double MinExtractedPhotonEnergy = 30.0; //MeV
  double FirstBurstTime  =      1000; //1e4
  double AverageInterval = 2000;//86400.0; //s
  bool  GenerateFluence  =   false;
  double Fluence = 1.0e-5;
  double Duration;
  int BURSTtype = 0; // 1->S, 2->L, else S+L
  bool GLASTCoordinate            = false;
  double theta = 45.0;
  double phi   = 0.0;

  int  GBM                        = 1;

  std::ofstream os("GRBobs_user_library.xml",std::ios::out);
  std::ofstream osTest("../src/test/GRBParam.txt",std::ios::out);

  os<<"<!-- $Header -->"<<std::endl;
  os<<"<!-- ************************************************************************** -->"<<std::endl;
  os<<"<source_library title=\"GRBobs_user_library\">"<<std::endl;
  char SourceName[100]="GRB";
  
  double BurstTime=FirstBurstTime;
  
  double alpha;
  double beta;
  osTest<<"    Seed Fluence Npeak  LowEnergy HighEnergy"<<std::endl;
  int type;
  int Nlong=0;
  int Nshort=0;
  for(int i = 0; i<N ; i++)
    {
      if(BURSTtype==1) type=1;
      else if(BURSTtype==2) type=2;
      else type = ((rnd->Uniform()<=0.3) ? 1 : 2);
      if (type==1)
	{
	  Duration = pow(10.0,(double)rnd->Gaus(-0.2,0.55)); //(Short Bursts)
	  Fluence = ((GenerateFluence) ? pow(10.0,(double)rnd->Gaus(-6.3,0.57)) : Fluence); //erg/cm^2 (Short Bursts)
	  Nshort++;
	}
      else
	{
	  Duration = pow(10.0,(double)rnd->Gaus(1.46,0.49)); //(Long Burst)
	  Fluence = ((GenerateFluence) ? pow(10.0,(double)rnd->Gaus(-5.4,0.62)): Fluence); //erg/cm^2 (Long Burst)
	  Nlong++;
	}

      theta = 0.0;//10.0*(i%8);
      phi   = 0.0;//10.0*(i%37);
      alpha= -1.0;//*rnd->Uniform(0.0, 1.5);
      beta = -2.25;//-1.0*rnd->Uniform(2.0, 2.5);
      os<<""<<std::endl;
      if(i<10) os<<"<source name=\" GRB_0000"<<i<<" \">"<<std::endl;
      else if(i<100) os<<"<source name=\" GRB_000"<<i<<" \">"<<std::endl;
      else if(i<1000) os<<"<source name=\" GRB_00"<<i<<" \">"<<std::endl;
      else if(i<10000) os<<"<source name=\" GRB_0"<<i<<" \">"<<std::endl;
      else os<<"<source name=\" GRB_"<<i<<" \">"<<std::endl;
      os<<"<spectrum escale=\"MeV\">"<<std::endl;
      os<<" <SpectrumClass name=\"GRBobsmanager\" params=\""<<BurstTime<<" , "<<Duration<<" , "<<Fluence<<" , "<<alpha<<" , "<<beta<<" , "<<MinExtractedPhotonEnergy<<" , "<<GBM<<" \"/>"<<std::endl;
      if(GLASTCoordinate)  
	os<<"<direction theta=\""<<theta<<"\" phi=\""<<phi<<"\" />"<<std::endl;//u galactic_dir l=\""<<l<<"\" b=\""<<b<<"\" />"<<std::endl;
      else
	os<<"<use_spectrum frame=\"galaxy\" />"<<std::endl;//u galactic_dir l=\""<<l<<"\" b=\""<<b<<"\" />"<<std::endl;
      os<<"</spectrum> </source>"<<std::endl;
      //////////////////////////////////////////////////
      osTest<<BurstTime<<" "<<Duration<<" "<<Fluence<<" "<<alpha<<" "<<beta<<std::endl;
      //////////////////////////////////////////////////
      BurstTime+=rnd->Exp(AverageInterval);
    }
  os<<""<<std::endl;
  os<<"<source name=\" AllGRBobs\" >"<<std::endl;
  for(int i = 0; i<N ; i++)
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



