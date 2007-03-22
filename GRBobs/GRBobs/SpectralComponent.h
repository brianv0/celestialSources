#ifndef SpectralComponent_h
#define SpectralComponent_h 1

#include "SpectObj/SpectObj.h"

class SpectralComponent
{
 public:
  SpectralComponent(SpectObj *spectrum,double tstart,double tend);
  
  double flux(double time, double MinPhotonEnergy);
  double interval(double time, double MinPhotonEnergy);
  double energyMeV(double time, double MinPhotonEnergy);
  inline void SetResiduals(double res) {m_residuals = res;}

 private:
  double m_tstart;
  double m_tend;
  double m_residuals;
  
  SpectObj *m_spectrum;
  
};
#endif


