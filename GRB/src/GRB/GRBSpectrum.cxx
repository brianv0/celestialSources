///
///   GRBSpectrum: Spectrum class for the GRB source
///    Authors: Nicola Omodei & Johann Cohen Tanugi 
///
#include "GRBSpectrum.h"
#include "CLHEP/Random/RandFlat.h"
#include <iostream>
#include <cmath>

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
	  std::cout<<"Restarting GRBSim..."<<std::endl;
	  flag=0;
	}
    }
}


///return flux, given a time
double GRBSpectrum::flux(double time) const
{
  return m_spectrum.integrated_Flux(m_grbsim->EnergyPh());//IRate(m_spectrum,m_grbsim->EnergyPh()); // in ph/(m^2 s) 
}

double GRBSpectrum::interval(double time)// const
{
  double inter;
  double tmax=m_grbsim->Tmax();
  double sum=0.0;         // Initialize the integral=0
  double temp;
  
  if(time<=1.0e-9) time=1.0e-9;
  temp=flux(time); // The m_spectrum is already computed from energySrc()
  //std::cout<<" Interval @ time "<<time<<" Rate = "<<temp<<std::endl;
  
    /// test to implement the right rate...
  double t1=time;         // t1 is the integration variable...
  double dt=1.0e-3*tmax;       // .. and dt is its integration step.
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
	temp=(m_grbsim->ComputeFlux(t1)).integrated_Flux(m_grbsim->EnergyPh()); // is the rate at t1...
	inter = t1-time;
      }
  }//end of while
  //  if (inter <= 1.0e-8) inter = 1.0e-8;
  return inter; // in Seconds 
}

double GRBSpectrum::energy(double time)
{
  std::cout<<"energySrc @ time "<<time<<std::endl;
  //return the flux in photons/(m^2 sec)
  /// Use time to update the spectrum
  //  m_spectrum.clear();
  //  m_grbsim->ComputeFlux(time);
  m_spectrum = m_grbsim->ComputeFlux(time);//;m_grbsim->Spectrum();
  //m_grbsim->setSpectrum(m_spectrum);
    /// Then returns a uniform random value, to be used by operator()
  return (*this)(RandFlat::shoot(1.0));
}

float GRBSpectrum::operator() (float u) const
{
  double energy = m_spectrum.DrawPhotonFromSpectrum(u,m_grbsim->EnergyPh());
  return (float) energy;
}




