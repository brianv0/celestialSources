#include "RadiationProcess.h"
#include "GRBConstants.h"

RadiationProcess::RadiationProcess(SpectObj spectrumObj)
  : m_spectrumObj(spectrumObj)
{;}

double RadiationProcess::processFlux(double E, 
				     double ec, 
				     double em)
{
  // The slopes in the different energies are taken from the synchrotron 
  // radiation.
  double value = 0;
  if(ec <= em) //FAST COOLING
    {
      if(E <= ec)
	value = pow(E/ec,0.33333);
      else if (E <= em)
	value = pow(E/ec,-0.5);
      else 
	value = pow(em/ec,-0.5)*pow(E/em,-cst::p/2.);
    }
  else //SLOW COOLING
    {
      if(E <= em)
	value = pow(E/em,0.33333);
      else if (E <= ec)
	value = pow(E/em,-(cst::p-1.)/2.);
      else 
	value = pow(ec/em,-(cst::p-1.)/2.)*pow(E/ec,-cst::p/2.);
    }
  return value;
}

double RadiationProcess::electronNumber(double gi,
					double gamma_min, 
					double gamma_max,
					double dr,
					double ComovingTime, 
					double CoolingTime,
					double N0)
{
  //if (gi < gamma_min) gi = gamma_min;
  
  //SPECTRAL DISTRIBUTION
  const double N_0  = N0*(cst::p-1.)/
    (pow(gamma_min,1.0-cst::p)-pow(gamma_max,1.0-cst::p)); 
  double N_e = N_0;
  
  if (gi < gamma_min || gi > gamma_max) 
    N_e *= 0.0;
  else 
    N_e *= pow(gi,-cst::p); //adim;
  
  // TEMPORAL DISTRIBUTION
  
  double gi2     = pow(gi,2.);
  
  if(cst::pulse_shape=="sgauss")
    {
      double tpeak   = (dr/cst::c)*sqrt((gi+1.)/(gi-1.));
      // Simmetric Pulse shape. 
      //The intrinsic delay is aprrox 0. Rise time = decay time 
      if (ComovingTime > 2.*tpeak)
 	N_e=0.;
      else
	N_e *= (exp(-abs(ComovingTime-tpeak)/CoolingTime)
		-exp(-tpeak/CoolingTime))/(1.-exp(-tpeak/CoolingTime));
    }
  else if(cst::pulse_shape=="agauss")
    {
      // Asimmetric Pulse shape. The lag between different channel is visible
      // Rise time != decay time. 
      //double tcross  = (dr/cst::c)*gamma_min/gi;
      double tpeak   = (dr/cst::c)*gamma_min/gi;
      
      if(ComovingTime < tpeak)
	N_e *= (exp((ComovingTime-tpeak)/tpeak)-exp(-1.))/(1.-exp(-1.));
      else
	N_e *= (exp(-(ComovingTime-tpeak)/CoolingTime)); 
      
      //        if(ComovingTime < tpeak)
      //  	N_e *= (exp((ComovingTime-tpeak)/tcross)-exp(-tpeak/tcross))/
      //        (1.-exp(-tpeak/tcross));
      //        else
      //  	N_e *= (exp(-(ComovingTime-tpeak)/CoolingTime)); 
    }
  else
    {
      // Exponential shape, FRED like function...
      
      double tpeak   = (dr/cst::c)/gamma_min*(gamma_min*gamma_min-1.)
	*gi/(gi2-1.);
      if(ComovingTime < tpeak)
	N_e *= CoolingTime/tpeak*(1.-exp(-ComovingTime/CoolingTime));
      else
	N_e *= CoolingTime/tpeak*(1.-exp(-tpeak/CoolingTime))*
	  exp((tpeak-ComovingTime)/CoolingTime);
    }
  return N_e;
} 

double RadiationProcess::timeShiftForDispersion(const double time, 
						const double E, 
						const double distance_to_source)
{
  // R in eV; converted in GeV
  const double E_QG = 1.e+28; //eV
  const double ksi = 1.;
  double dispersion_factor = ksi * E/ E_QG;
  double shifted_time = time-(distance_to_source / cst::c)*dispersion_factor;
  return shifted_time;
}


double RadiationProcess::comovingTime(const double time, 
				      const double gamma, 
				      const double E, 
				      const double distance_to_source)
{
  if(cst::flagQG)
    return gamma * timeShiftForDispersion(time,E, distance_to_source);
  else 
    return gamma * time;
}
