//////////////////////////////////////////////////
void GenerateXMLLibrary(int N=100)
{
  TRandom *rnd = new TRandom();
  double MinExtractedPhotonEnergy = 30.0; //MeV
  double FirstBurstTime  =      10000;
  double AverageInterval = 86400.0; //s
  bool  GenerateFluence  =   true;
  double Fluence = 1.0e-5;
  
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

  int Npulses;
  double alpha;
  double beta;
  osTest<<"    Seed Fluence Npeak  LowEnergy HighEnergy"<<std::endl;

  for(int i = 0; i<N ; i++)
    {
      int Npulses = rnd->Uniform(1,10);
      theta = 0.0;//10.0*(i%8);
      phi   = 0.0;//10.0*(i%37);
      alpha= -1.0*rnd->Uniform(0.0, 1.5);
      beta = -1.0*rnd->Uniform(2.0, 2.5);
      // L Fluence = pow(10.0,rnd->Gaus(-5.4,0.63));
      // S Fluence = pow(10.0,rnd->Gaus(-6.4,0.57));
      Fluence = (GenerateFluence) ? pow(10.0,rnd->Gaus(-5.6,0.76)) : Fluence;
      os<<""<<std::endl;
      if(i<10) os<<"<source name=\" GRB_000"<<i<<" \">"<<std::endl;
      else if(i<100) os<<"<source name=\" GRB_00"<<i<<" \">"<<std::endl;
      else if(i<1000) os<<"<source name=\" GRB_0"<<i<<" \">"<<std::endl;
      else if(i<1000) os<<"<source name=\" GRB_"<<i<<" \">"<<std::endl;
      else os<<"<source name=\" GRB_0"<<i<<" \">"<<std::endl;
      os<<"<spectrum escale=\"MeV\">"<<std::endl;
      os<<" <SpectrumClass name=\"GRBobsmanager\" params=\""<<BurstTime<<" , "<<Fluence<<" , "<<Npulses<<" , "<<alpha<<" , "<<beta<<" , "<<MinExtractedPhotonEnergy<<" , "<<GBM<<" \"/>"<<std::endl;
      if(GLASTCoordinate)  
	os<<"<direction theta=\""<<theta<<"\" phi=\""<<phi<<"\" />"<<std::endl;//u galactic_dir l=\""<<l<<"\" b=\""<<b<<"\" />"<<std::endl;
      else
	os<<"<use_spectrum frame=\"galaxy\" />"<<std::endl;//u galactic_dir l=\""<<l<<"\" b=\""<<b<<"\" />"<<std::endl;
      os<<"</spectrum> </source>"<<std::endl;
      //////////////////////////////////////////////////
      osTest<<BurstTime<<" "<<Fluence<<" "<<Npulses<<" "<<alpha<<" "<<beta<<std::endl;
      //////////////////////////////////////////////////
      BurstTime+=rnd->Exp(AverageInterval);
    }
  os<<""<<std::endl;
  os<<"<source name=\" AllGRBobs\" >"<<std::endl;
  for(int i = 0; i<N ; i++)
    {
      if(i<10) os<<"     <nestedSource sourceRef=\" GRB_000"<<i<<" \"/>"<<std::endl;
      else if(i<100) os<<"     <nestedSource sourceRef=\" GRB_00"<<i<<" \"/>"<<std::endl;
      else if(i<1000) os<<"     <nestedSource sourceRef=\" GRB_0"<<i<<" \"/>"<<std::endl;
      else if(i<1000) os<<"     <nestedSource sourceRef=\" GRB_0"<<i<<" \"/>"<<std::endl;
      else os<<"     <nestedSource sourceRef=\" GRB_"<<i<<" \"/>"<<std::endl;
    }
  os<<"</source>"<<std::endl;
  os<<"  </source_library>"<<std::endl;
}



