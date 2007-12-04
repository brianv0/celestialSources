void GenerateFakeBurst(int NumberOfPhotonsAbove100MeV = 10) // per square meter
{
  using namespace std;
  int EnergyBins = 500;
  double MinimumEnergy = 10.0;
  double MaximumEnergy = 3.0e8;

  int TimeBins    = 1000;
  double DeltaT      = 0.016;

  double de   = pow(MaximumEnergy/MinimumEnergy,1.0/(EnergyBins-1)); 

  ofstream of("test.dat");
  of<<"EnergyBins= "<<EnergyBins<<endl;
  of<<"MinimumEnergy= "<<MinimumEnergy<<endl;
  of<<"MaximumEnergy= "<<MaximumEnergy<<endl;
  of<<"TimeBins= "<<TimeBins<<endl;
  of<<"TimeBinWidth= "<<DeltaT<<endl;
  

  double energy, time, Nv, pulse;
  double t0 = TimeBins*DeltaT/2.;
  double a = 2.25;  
  //  double e1 = 3.0e4;
  double e1 = 1.0e5;
  double e2 = MaximumEnergy;


  //  double Const = 1e4*2.0*pow(t0,2.)*(e2-e1)*(pow(e1,1.-a) - pow(e2,1.-a))/(e1-a*e1-e2+a*e2);
  double Const = pow(t0,2.)*(pow(e1,1.-a) - pow(e2,1.-a))/(a-1.);
  for(int ti=0; ti<TimeBins;ti++)
    {
      time = ti*DeltaT;
      if(time<t0) 
	pulse = time;
      else 
	pulse = 2.*t0-time;
      
      for(int ei=0; ei<EnergyBins;ei++)
	{
	  energy = MinimumEnergy * pow(de,1.0*ei); //keV
	  
	  Nv = (NumberOfPhotonsAbove100MeV)/(1e4*Const)*TMath::Max(0.0,pulse)*pow(energy,-a);// [ph/(cm² s keV)]
	  of<<Nv<<"\t";
	}
      of<<endl;
    }
  of.close();
}
