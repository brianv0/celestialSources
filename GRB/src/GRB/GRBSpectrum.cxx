///
///   GRBSpectrum: Spectrum class for the GRB source
///    Authors: Nicola Omodei & Johann Cohen Tanugi 
///
#include "GRBSpectrum.h"
#include "CLHEP/Random/RandFlat.h"
#include <iostream>
#include <math.h>

GRBSpectrum::GRBSpectrum(const std::string& params) 
{
  m_grbsim = new GRBSim();
  int flag=0;
  while(flag==0)
    {
      try 
	{ 
	  m_grbsim->Start();
	  flag=1;
	}
      catch (char * s)
	{
	  std::cout<<"Restarting GRBSim...\n";
	  flag=0;
	}
    }
}

GRBSpectrum::~GRBSpectrum() 
{
  delete m_grbsim;
}

// Must be set to 1.0 for a point source.
double GRBSpectrum::solidAngle() const
{
  return 1.0;
}

///return flux, given a time
double GRBSpectrum::flux(double time) const
{
  //  double rateout;
  m_grbsim->ComputeFlux(time);
  /// test to implement the right rate...
  // if (m_grbsim->IRate()<=0.1?rateout=0.1:rateout=m_grbsim->IRate());
  //  cout<< " time ="<<time<<" Rate = "<<m_grbsim->IRate()<< endl;
  return m_grbsim->IRate(); // in ph/(m^2 s) 
}

double GRBSpectrum::rate(double time) const
{
  m_grbsim->ComputeFlux(time);
  /// test to implement the right rate...
  //  cout<< " time ="<<time<<" Rate = "<<m_grbsim->IRate()<< endl;
  return m_grbsim->IRate(); // in ph/(m^2 s) 
}

double GRBSpectrum::interval(double time)// const
{
  //cout<<" INTERVAL!<<endl;
  /// test to implement the right rate...
  double sum=0.0;
  double temp;
  double t1=time;
  double dt;
  
  while (sum<1.0){
    dt=1.0e-2;
    temp=rate(t1);
    if (temp>1/cst::DeadTime){return t1-time+cst::DeadTime;}
    if (temp<=1.0) {temp=1.0;}
    else if (temp>1.0/dt) {return t1-time+1/temp;}
    sum+=temp*dt;
    t1+=dt;
  }
  return t1-time; // in Seconds 
}

double GRBSpectrum::energySrc(HepRandomEngine* engine, double time)
{
  //return the flux in photons/(m^2 sec)
  /// Use time to update the spectrum
  m_spectrum.clear();
  m_grbsim->ComputeFlux(time);
  m_spectrum =m_grbsim->Spectrum();
  /// Then returns a uniform random value, to be used by operator()
  return (*this)(engine->flat());
}

float GRBSpectrum::operator() (float u) const
{
 
  float energy = m_grbsim->DrawPhotonFromSpectrum(m_spectrum, u);
  return (float) energy;
}

std::pair<float,float> GRBSpectrum::dir(float energy) const
{  
  return m_grbsim->GRBdir();
}



