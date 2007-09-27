#include <iomanip>
//////////////////////////////////////////////////
// OPTIONS FOR GENERATING THE XML LIBRARY:
double MinExtractedPhotonEnergy = 10.0; //MeV
long FirstBurstTime    =      1000; //1e4
double AverageInterval =      1000;
TString GenerateIC       =   "no";//"random"; //"yes", "no"
bool  GenerateGBM        =   true;//false;
bool  GeneratePF        =   false;

//////////////////////////////////////////////////

double alpha     = -1.001;
double Epeak     = 500.0;
double phi       = 0.0;
double z=0;

//////////////////////////////////////////////////

double t90_min   = 1.0;
double t90_max   = 100.0;
int    NlogT90   = 3;


double beta_min  = -3.2;
double beta_max  = -1.2;
int    NBeta     = 5;
double deltaBeta = (beta_max-beta_min)/(NBeta-1);

double PF_min    =  0.1;
double PF_max    = 30.0;
int    NlogPF       = 5;

double FL_min    =  1e-8;
double FL_max    =  1e-4;
int    NlogFL    =  8;

double theta_min  = 0.0;
double theta_max  = 90.0;
int    Ntheta     = 5;
double Dtheta     = (theta_max - theta_min)/(Ntheta-1);
//////////////////////////////////////////////////
long seed = 1;
double Essc_Esyn=1.0e4;
double Fssc_Fsyn=0.0;
double co_energy=0.0;

void GenerateGRBGrid()
{
  int Nbursts = 0;
  if(GeneratePF)
    Nbursts=NlogT90*NBeta*NlogPF*Ntheta;
  else
    Nbursts=NlogT90*NBeta*NlogFL*Ntheta;
  
  std::cout<<" Number of bursts: "<<Nbursts<<std::endl;
  
  //////////////////////////////////////////////////
  std::ofstream osXML("GRBobs_grid_library.xml",std::ios::out);
  std::ofstream osTest("../src/test/GRBParam.txt",std::ios::out);
  osXML.precision(4);
  
  osXML<<"<!-- $Header -->"<<std::endl;
  osXML<<"<!-- ************************************************************************** -->"<<std::endl;
  osXML<<"<source_library title=\"GRBobs_user_library\">"<<std::endl;
  char SourceName[100]="GRB";
  long BurstTime=FirstBurstTime;
  
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
  int NGRBs=0;
  for(int t = 0; t<Ntheta ; t++)
    {
      double theta = theta_min + t * Dtheta;
      for(int j = 0; j<NlogT90; j++)
	{
	  double Duration = t90_min * pow(t90_max/t90_min ,1.0*j/(NlogT90-1));
	  for(int l = 0; l<NBeta; l++)
	    {
	      double beta = beta_min + deltaBeta * l;
	      if(GeneratePF)
		{
		  for(int m = 0; m<NlogPF; m++)
		    {
		      double PeakFlux = PF_min * pow(PF_max/PF_min,1.0*m/(NlogPF-1));
		      osXML<<""<<std::endl;
		      if(NGRBs<10) osXML<<"<source name=\" GRB_0000"<<NGRBs<<" \">"<<std::endl;
		      else if(NGRBs<100) osXML<<"<source name=\" GRB_000"<<NGRBs<<" \">"<<std::endl;
		      else if(NGRBs<1000) osXML<<"<source name=\" GRB_00"<<NGRBs<<" \">"<<std::endl;
		      else if(NGRBs<10000) osXML<<"<source name=\" GRB_0"<<NGRBs<<" \">"<<std::endl;
		      else osXML<<"<source name=\" GRB_"<<NGRBs<<" \">"<<std::endl;
		      osXML<<"<spectrum escale=\"MeV\">"<<std::endl;
		      osXML<<" <SpectrumClass name=\"GRBobsmanager\" params=\"seed="<<seed;
		      osXML<<", tstart="<<BurstTime<<", duration="<<Duration<<", peakFlux="<<PeakFlux;
		      osXML<<", theta="<<theta<<", phi="<<phi;
		      osXML<<", redshift="<<z<<", alpha="<<alpha<<", beta="<<beta<<", Ep="<<Epeak<<", emin="<<MinExtractedPhotonEnergy;
		      osXML<<", GBM="<<(int)GenerateGBM;
		      osXML<<" \"/>"<<std::endl;
		      osXML<<"<use_spectrum frame=\"galaxy\" />"<<std::endl; //u galactic_dir l=\""<<l<<"\" b=\""<<b<<"\" />"<<std::endl;
		      osXML<<"</spectrum> </source>"<<std::endl;
		      osTest<<setw(11)<<BurstTime<<" "<<setw(8)<<setprecision(4)<<Duration<<" "<<setw(8)<<PeakFlux<<" ";
		      osTest<<setw(8)<<z<<" "<<setw(8)<<alpha<<" "<<setw(8)<<beta<<" "<<setw(8)<<Epeak<<" ";
		      osTest<<setw(8)<<Essc_Esyn<<setw(8)<<Fssc_Fsyn<<setw(8)<<co_energy<<std::endl;
		      //////////////////////////////////////////////////
		      BurstTime+=(int) (AverageInterval);		    
		      NGRBs++;		  
		    }
		}
	      else
		{ 
		  for(int m = 0; m<NlogFL; m++)
		    {
		      double fluence = FL_min * pow(FL_max/FL_min,1.0*m/(NlogFL-1));
		      osXML<<""<<std::endl;
		      if(NGRBs<10) osXML<<"<source name=\" GRB_0000"<<NGRBs<<" \">"<<std::endl;
		      else if(NGRBs<100) osXML<<"<source name=\" GRB_000"<<NGRBs<<" \">"<<std::endl;
		      else if(NGRBs<1000) osXML<<"<source name=\" GRB_00"<<NGRBs<<" \">"<<std::endl;
		      else if(NGRBs<10000) osXML<<"<source name=\" GRB_0"<<NGRBs<<" \">"<<std::endl;
		      else osXML<<"<source name=\" GRB_"<<NGRBs<<" \">"<<std::endl;
		      osXML<<"<spectrum escale=\"MeV\">"<<std::endl;
		      osXML<<" <SpectrumClass name=\"GRBobsmanager\" params=\"seed="<<seed;
		      osXML<<", tstart="<<BurstTime<<", duration="<<Duration<<", fluence="<<fluence;
		      osXML<<", theta="<<theta<<", phi="<<phi;
		      osXML<<", redshift="<<z<<", alpha="<<alpha<<", beta="<<beta<<", Ep="<<Epeak<<", emin="<<MinExtractedPhotonEnergy;
		      osXML<<", GBM="<<(int)GenerateGBM;
		      osXML<<" \"/>"<<std::endl;
		      osXML<<"<use_spectrum frame=\"galaxy\" />"<<std::endl; //u galactic_dir l=\""<<l<<"\" b=\""<<b<<"\" />"<<std::endl;
		      osXML<<"</spectrum> </source>"<<std::endl;
		      osTest<<setw(11)<<BurstTime<<" "<<setw(8)<<setprecision(4)<<Duration<<" "<<setw(8)<<fluence<<" ";
		      osTest<<setw(8)<<z<<" "<<setw(8)<<alpha<<" "<<setw(8)<<beta<<" "<<setw(8)<<Epeak<<" ";
		      osTest<<setw(8)<<Essc_Esyn<<setw(8)<<Fssc_Fsyn<<setw(8)<<co_energy<<std::endl;
		      
		      //////////////////////////////////////////////////
		      BurstTime+=(int) (AverageInterval);		    
		      NGRBs++;
		    }
		}
	    }
	}
    }
  osXML<<""<<std::endl;
  osXML<<"<source name=\" AllGRBobs\" >"<<std::endl;
  for(int i = 0; i<NGRBs ; i++)
    {
      if(i<10) osXML<<"     <nestedSource sourceRef=\" GRB_0000"<<i<<" \"/>"<<std::endl;
      else if(i<100) osXML<<"     <nestedSource sourceRef=\" GRB_000"<<i<<" \"/>"<<std::endl;
      else if(i<1000) osXML<<"     <nestedSource sourceRef=\" GRB_00"<<i<<" \"/>"<<std::endl;
      else if(i<10000) osXML<<"     <nestedSource sourceRef=\" GRB_0"<<i<<" \"/>"<<std::endl;
      else osXML<<"     <nestedSource sourceRef=\" GRB_"<<i<<" \"/>"<<std::endl;
    }
  osXML<<"</source>"<<std::endl;
  osXML<<"  </source_library>"<<std::endl;
  
  std::cout<<"--------------------------------------------------"<<std::endl;
  std::cout<<"Generate "<<NGRBs<<" GRBs"<<std::endl;
  std::cout<<"--------------------------------------------------"<<std::endl;
  //  if (GenerateRedshift)  Redshifts->Print("GeneratedRedshift.eps");


}



