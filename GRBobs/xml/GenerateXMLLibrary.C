//////////////////////////////////////////////////
void GenerateXMLLibrary(int N=100)
{
  TRandom *rnd = new TRandom();
  double MinExtractedPhotonEnergy = 30.0; //MeV
  double FirstBurstTime=0;
  double AverageInterval=1000.0; //s
  
  bool GLASTCoordinate            = false;  
  std::ofstream os("GRBobs_user_library.xml",std::ios::out);
  std::ofstream osTest("../src/test/GRBobsParam.txt",std::ios::out);

  os<<"<!-- $Header -->"<<std::endl;
  os<<"<!-- ************************************************************************** -->"<<std::endl;
  os<<"<source_library title=\"GRBobs_user_library\">"<<std::endl;
  char SourceName[100]="GRB";
  
  double BurstTime=FirstBurstTime;
  double Fluence;
  int Npulses;
  double alpha;
  double beta;
  osTest<<"    Seed l     b   Fluence Npeak  LowEnergy HighEnergy"<<std::endl;

  
  for(int i = 0; i<N ; i++)
    {
      double l = rnd->Uniform(-180.0,180.0);
      double b = rnd->Uniform(-90.0,90.0);
      int Npulses = rnd->Uniform(1,10);
      alpha= -1.0;
      beta = -2.0;
      Fluence = pow(10.0,rnd->Gaus(-5.4,0.62));
      os<<""<<std::endl;
      os<<"<source name=\" GRB_"<<i<<" \">"<<std::endl;
      os<<"<spectrum escale=\"MeV\">"<<std::endl;
      os<<" <SpectrumClass name=\"GRBobsmanager\" params=\""<<l<<" , "<<b<<" ,"<<BurstTime<<" , "<<Fluence<<" , "<<Npulses<<" , "<<alpha<<" , "<<beta<<" , "<<MinExtractedPhotonEnergy<<"\"/>"<<std::endl;
      if(GLASTCoordinate)  
	os<<"<direction theta=\"45\" phi=\"0\" />"<<std::endl;//u galactic_dir l=\""<<l<<"\" b=\""<<b<<"\" />"<<std::endl;
      else
	os<<"<use_spectrum frame=\"galaxy\" />"<<std::endl;//u galactic_dir l=\""<<l<<"\" b=\""<<b<<"\" />"<<std::endl;
      os<<"</spectrum> </source>"<<std::endl;
      //////////////////////////////////////////////////
      osTest<<Npulses<<" "<<l<<" "<<b<<" "<<Fluence<<" "<<Npulses<<" "<<alpha<<" "<<beta<<std::endl;
      //////////////////////////////////////////////////
      BurstTime+=rnd->Exp(AverageInterval);
    }
  os<<""<<std::endl;
  os<<"<source name=\" AllGRBobs\" >"<<std::endl;
  for(int i = 0; i<N ; i++)
    {
      os<<"     <nestedSource sourceRef=\" GRB_"<<i<<"\" />"<<std::endl;
    }
  os<<"</source>"<<std::endl;
  os<<"  </source_library>"<<std::endl;
}



