
#ifndef RadiationProcess_H
#define RadiationProcess_H 1

#include "SpectObj.h"
#include"GRBConstants.h"

class RadiationProcess
{

 public:
  RadiationProcess(){;}
  RadiationProcess(SpectObj);
  virtual ~RadiationProcess(){;}

  virtual double processFlux(double, double, double);
  virtual double electronDensity(double, double, double,
				 double, double, double, 
				 double, double);

  //! return the Spectrum object.
  inline SpectObj getSpectrumObj() {return m_spectrumObj;}

  //! Helper function. /param B is in Gauss, Umag is in erg.
  inline virtual double Umag(const double B)
    {
      return pow(B/1.0e+4,2.)/(200.0*M_PI*cst::mu0);
    }
 
 protected:
  SpectObj m_spectrumObj;

};

#endif;
