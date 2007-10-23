#include "GRBobs/SpectralComponent.h"

SpectralComponent::SpectralComponent(SpectObj *spectrum,double tstart,double tend)
{
  m_tstart = tstart;
  m_tend   = tend;
  m_spectrum = spectrum;
  m_residuals = 0.0;
}

double SpectralComponent::flux(double time, double MinPhotonEnergy)
{
  double flux=0.0;
  if (time >= m_tstart && time <= m_tend) 
    flux = m_spectrum->flux(time-m_tstart,MinPhotonEnergy);
  //  std::cout<<"flux: "<<flux<<" t = "<<time<<" ts = "<<m_tstart<<" te = "<<m_tend<<std::endl;
  return flux;  
}

double SpectralComponent::interval(double time, double MinPhotonEnergy)
{
  if (time<m_tstart) 
    return m_tstart - time;
  if(time>m_tend) 
    return 1e10;
  
  double inte = m_spectrum->interval(time-m_tstart,MinPhotonEnergy) - m_residuals;
  //  std::cout<<"inte: "<<inte<<" ts = "<<m_tstart<<" te = "<<m_tend<<std::endl;
  return inte;
}

double SpectralComponent::energyMeV(double time, double MinPhotonEnergy)
{
  double ene=0.1;
  if(m_residuals==0)
    {
      ene =  m_spectrum->energy(time-m_tstart,MinPhotonEnergy)*1.0e-3; //MeV
      //      std::cout<<"ene( "<<time<<")= "<<ene<<std::endl;
    }
  return ene;
}
