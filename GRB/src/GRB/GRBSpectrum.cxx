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
  m_grbsim = new GRBSim(params);
  int flag=0;
  while(flag==0)
    {
      try 
	{ 
	  m_grbsim->MakeGRB();
	  flag=1;
	}
      catch (char * s)
	{
	  std::cout<<"Restarting GRBSim...\n";
	  flag=0;
	}
    }
}


///return flux, given a time
double GRBSpectrum::flux(double time) const
{
  //  cout<<"FLUX"<<endl;
  //m_grbsim->ComputeFlux(time);
  /// test to implement the right rate...
  return m_grbsim->IRate(m_spectrum,m_grbsim->EnergyPh()); // in ph/(m^2 s) 
}

double GRBSpectrum::interval(double time)// const
{
  //  cout<<" Interval @ time "<<time<<endl;;
  double inter;
  if(time<0) time=abs(time);
  
  double tmax=m_grbsim->Tmax();
  /// test to implement the right rate...
  double sum=0.0;         // Initialize the integral=0
  double t1=time;         // t1 is the integration variable...
  double temp=flux(time); // The m_spectrum is already computed from energySrc()
  double dt=1.0e-2;       // .. and dt is its integration step.
  while (sum<1.0){ // It compute the integral of the rate...
    // The idea is to compute t1 for which the 
    // integral between time and t1 of the rate(t) dt =1.
    if(t1>=tmax)
      {
	inter = 1.0e+6;
	break; // This exits the loop in the case af the burst is finished
      } 
    if (temp<=1/tmax) temp=1/tmax; // Minumum rate...
    if (temp>=1.0/dt)  // the sum will be > 1, and I can stop the while loop...
      {
	inter = t1-time+1.0/temp;
	break;
      }
    else
      {
	sum+=temp*dt; //...integrate...
	t1+=dt;       // ... and  incrase the time...
	temp=m_grbsim->IRate(m_grbsim->ComputeFlux(t1),
			     m_grbsim->EnergyPh()); // is the rate at t1...
	inter = t1-time;
      }
  }//end of while
  //  if (inter <= 1.0e-8) inter = 1.0e-8;
  return inter; // in Seconds 
}

double GRBSpectrum::energySrc(HepRandomEngine* engine, double time)
{
  //  cout<<"energySrc @ time "<<time<<endl;
  //return the flux in photons/(m^2 sec)
  /// Use time to update the spectrum
  m_spectrum.clear();
  //  m_grbsim->ComputeFlux(time);
  m_spectrum = m_grbsim->ComputeFlux(time);//;m_grbsim->Spectrum();
  m_grbsim->setSpectrum( m_spectrum);
    /// Then returns a uniform random value, to be used by operator()
  return (*this)(engine->flat());
}

float GRBSpectrum::operator() (float u) const
{
  
  double energy = m_grbsim->DrawPhotonFromSpectrum(m_spectrum, 
						   u,m_grbsim->EnergyPh());
  return (float) energy;
}

std::pair<float,float> GRBSpectrum::dir(float energy) const
{  
  return m_grbsim->GRBdir();
}



