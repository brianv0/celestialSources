/**
* @file EarthPhenomLimb.cxx
* @brief A phenomenological model of the Earth limb based on LAT measurements
* @author Keith Bechtol
*
* $Header$
*/

#include <iostream>

#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <vector>
#include <stdlib.h>

#include "facilities/Util.h"

#include "flux/SpectrumFactory.h"

#include "EarthPhenom/EarthPhenomLimb.h"
#include "TString.h"
#include "TGraph.h"

#include "CLHEP/Random/RandFlat.h"

#include "celestialSources/ConstParMap.h"

//static SpectrumFactory<EarthPhenomLimb> factory;
//const ISpectrumFactory& EarthPhenomLimbFactory = factory;

EarthPhenomLimb::EarthPhenomLimb(const std::string &paramString)
  : m_normalization(1.), m_emin(10.), m_emax(350000.), m_invert_direction(false) {	
  
  celestialSources::ConstParMap pars(paramString);
  
  m_normalization = pars.value("norm"); // Normalization = 1 is default
  m_emin = pars.value("emin"); // MeV
  m_emax = pars.value("emax"); // MeV
  if(pars.value("invert_zenith") == 1) // 0 = normal mode, 1 = invert direction
    m_invert_direction = true;
  else
    m_invert_direction = false;

  init(m_normalization, m_emin, m_emax);

  std::cerr << "EarthPhenomLimb (updated August 2013) created. Normalization = "
	    << m_normalization << " . Total flux = "
	    << m_integral_flux << " m^-2 s^-1" << " between "
	    << m_emin << " MeV and "
	    << m_emax << " MeV." << std::endl;

  if(m_invert_direction)
    std::cerr << "Running with inverted photon directions..." << std::endl;
  else
    std::cerr << "Running in normal mode..." << std::endl;

  m_energy_called = false;
}

double EarthPhenomLimb::flux(double time) const { // Argument is the mission elapsed time (s)
    return m_integral_flux / solidAngle(); // (m^-2 s^-1 sr^-1)
}

double EarthPhenomLimb::solidAngle() const {
  return m_solid_angle; // (sr)
}

void EarthPhenomLimb::init(double normalization, double emin, double emax) {
  
  // Initialize pre-fitted model parameters

  // Spectral component
  m_spectral_prefactor = normalization * 1.e4 * 1.031e-01; // (MeV^-1 m^-2 s^-1 sr^-1) Converting from cm^-2 to m^-2
  m_spectral_index1 = -1.532;
  m_spectral_index2 = -2.790;
  m_spectral_ebreak = 3.703e+02; // (MeV)
  m_spectral_beta = 7.276e-01;

  // Zenith angle component
  // Zenith = 0 deg -> away from Earth center; zenith = 180 deg -> towards Earth center
  m_zenithmin=180.-70.; // (deg)
  m_zenithmax=180.-50.; // (deg)
  m_zenith_peak = 180.-6.795e+01; // (deg) Converting from Earth nadir angle to Earth zenith angle 
  m_zenith_width = 3.184e-01; // (deg) Deconvolve using a 68% containment radius of 0.17 deg
  
  m_zenith_energy_slope_prefactor = 5.401e-04;
  m_zenith_energy_slope_index = 9.272e-01;
  m_zenith_energy_slope_mu = 1.000e+03; // (MeV)
  m_zenith_energy_slope_tau = 4.401e+03; // (MeV)

  m_zenith_energy_curve_tanh_prefactor = -3.839e+00;
  m_zenith_energy_curve_tanh_cross = 4.000e+00; // log10(MeV)
  m_zenith_energy_curve_tanh_width = 1.399e-01; // (dex)
  m_zenith_energy_curve_pol_4 = 3.867e-04;

  // Azimuth angle component
  // Looking towards the Earth
  // LAT analysis coordinate system: Azimuth = 0 deg -> north; azimuth = 90 deg -> east (left-handed)
  m_azimuthmin=0.; // (deg)
  m_azimuthmax=360.; // (deg)
  m_azimuthal_logsine_phase = 176.5; // (deg)
  m_azimuthal_notch_phase = 106.0; // (deg)

  m_azimuthal_energy_logsine_prefactor = 6.323e-02;
  m_azimuthal_energy_logsine_index = 3.742e-01;
  m_azimuthal_energy_logsine_ecutoff = 7.133e+03; // (MeV)

  m_azimuthal_energy_notch_ratio_prefactor = 3.628e-03;
  m_azimuthal_energy_notch_ratio_index = 8.960e-01;
  m_azimuthal_energy_notch_ratio_ecutoff = 2.796e+02; // (MeV)
  m_azimuthal_energy_notch_width_prefactor = 8.871e-01;
  m_azimuthal_energy_notch_width_index = 7.975e-01;
  m_azimuthal_energy_notch_width_ecutoff = 9.283e+02; // (MeV)

  // Initialize model formulae

  // Variable is energy (MeV, linear-space)
  TString spectralEnergyFormula = "[0]*x^[1]*(1+(x/[3])^(([1]-[2])/[4]))^(-1.*[4])";
  m_fSpectralEnergy = TF1("fSpectralEnergy", spectralEnergyFormula.Data(), m_emin, m_emax);
  m_fSpectralEnergy.SetParameter(0,m_spectral_prefactor);
  m_fSpectralEnergy.SetParameter(1,m_spectral_index1);
  m_fSpectralEnergy.SetParameter(2,m_spectral_index2);
  m_fSpectralEnergy.SetParameter(3,m_spectral_ebreak);
  m_fSpectralEnergy.SetParameter(4,m_spectral_beta);

  // Variable is energy (MeV, log10-space)
  TString spectralFormula = "[0]*(10^x)^[1]*(1+((10^x)/[3])^(([1]-[2])/[4]))^(-1.*[4])*(pow(10,x)*log(10))"; // Weighted for probability
  m_fSpectral = TF1("fSpectral", spectralFormula.Data(), log10(m_emin), log10(m_emax));
  m_fSpectral.SetParameter(0,m_spectral_prefactor);
  m_fSpectral.SetParameter(1,m_spectral_index1);
  m_fSpectral.SetParameter(2,m_spectral_index2);
  m_fSpectral.SetParameter(3,m_spectral_ebreak);
  m_fSpectral.SetParameter(4,m_spectral_beta);

  // Variable is energy (MeV, log10-space)
  // Energy-dependent variation in inner Earth logarithmic slope
  TString zenithSlopeFormula = "[0]*(10^x)^[1]*exp(-1*(10^x-[2])/[3])";
  m_fZenithSlope = TF1("fZenithSlope", zenithSlopeFormula.Data(), log10(m_emin), log10(m_emax));
  m_fZenithSlope.SetParameter(0,m_zenith_energy_slope_prefactor);
  m_fZenithSlope.SetParameter(1,m_zenith_energy_slope_index);
  m_fZenithSlope.SetParameter(2,m_zenith_energy_slope_mu);
  m_fZenithSlope.SetParameter(3,m_zenith_energy_slope_tau);

  // Variable is energy (MeV, log10-space)
  // Energy-dependent variation in inner Earth logarithmic curve
  TString zenithCurveFormula = "0.5*[0]*(1+TMath::TanH((x-[1])/[2]))+(x<4)*([3]*(x-1)^4)";
  m_fZenithCurve = TF1("fZenithCurve", zenithCurveFormula.Data(), log10(m_emin), log10(m_emax));
  m_fZenithCurve.SetParameter(0,m_zenith_energy_curve_tanh_prefactor);
  m_fZenithCurve.SetParameter(1,m_zenith_energy_curve_tanh_cross);
  m_fZenithCurve.SetParameter(2,m_zenith_energy_curve_tanh_width);
  m_fZenithCurve.SetParameter(3,m_zenith_energy_curve_pol_4);
  
  // Variable is zenith angle (deg)
  // Parameters 2 and 3 are energy-dependent
  // Notice that zenith angle is used instead of nadir angle (sign flip)
  TString zenithFormula = "(x>[0])*exp([2]*([0]-x)+[3]*([0]-x)^2)+(x<=[0])*exp(-1*([0]-x)^2/(2.*[1]^2))";
  m_fZenith = TF1("fZenith", zenithFormula.Data(), m_zenithmin, m_zenithmax);
  m_fZenith.SetParameter(0,m_zenith_peak);
  m_fZenith.SetParameter(1,m_zenith_width);

  // Variable is energy (MeV, log10-space)
  // Energy-dependent variation of logsine component amplitude
  TString azimuthalLogsineFormula = "[0]*(10^x)^[1]*exp(-1*(10^x)/[2])";
  m_fAzimuthalLogsine = TF1("fAzimuthalLogsine", azimuthalLogsineFormula.Data(), log10(m_emin), log10(m_emax));
  m_fAzimuthalLogsine.SetParameter(0,m_azimuthal_energy_logsine_prefactor);
  m_fAzimuthalLogsine.SetParameter(1,m_azimuthal_energy_logsine_index);
  m_fAzimuthalLogsine.SetParameter(2,m_azimuthal_energy_logsine_ecutoff);

  // Variable is energy (MeV, log10-space)
  // Energy-dependent vairation of notch component amplitude relative to log-size component
  TString azimuthalNotchRatioFormula = "[0]*(10^x)^[1]*exp(-1*(10^x)/[2])";
  m_fAzimuthalNotchRatio = TF1("fAzimuthalNotchRatio", azimuthalNotchRatioFormula.Data(), log10(m_emin), log10(m_emax));
  m_fAzimuthalNotchRatio.SetParameter(0,m_azimuthal_energy_notch_ratio_prefactor);
  m_fAzimuthalNotchRatio.SetParameter(1,m_azimuthal_energy_notch_ratio_index);
  m_fAzimuthalNotchRatio.SetParameter(2,m_azimuthal_energy_notch_ratio_ecutoff);

  // Variable is energy (MeV, log10-space)
  // Energy-dependent vairation of notch component width
  TString azimuthalNotchWidthFormula = "[0]*(10^x)^[1]*exp(-1*(10^x)/[2])";
  m_fAzimuthalNotchWidth = TF1("fAzimuthalNotchWidth", azimuthalNotchWidthFormula.Data(), log10(m_emin), log10(m_emax));
  m_fAzimuthalNotchWidth.SetParameter(0,m_azimuthal_energy_notch_width_prefactor);
  m_fAzimuthalNotchWidth.SetParameter(1,m_azimuthal_energy_notch_width_index);
  m_fAzimuthalNotchWidth.SetParameter(2,m_azimuthal_energy_notch_width_ecutoff);

  // Variable is azimuth angle (deg)
  // Parameters 2, 3, and 4 are energy-dependent
  TString azimuthalFormula = "exp([2]*sin(TMath::Pi()*(x-[0])/180))-[3]*exp(-1*(x-[1]-360)**2/(2*pow([4],2)))-[3]*exp(-1*(x-[1])**2/(2*pow([4],2)))-[3]*exp(-1*(x-[1]+360)**2/(2*pow([4],2)))";
  m_fAzimuth = TF1("fAzimuth", azimuthalFormula.Data(), m_azimuthmin, m_azimuthmax);
  m_fAzimuth.SetParameter(0,m_azimuthal_logsine_phase);
  m_fAzimuth.SetParameter(1,m_azimuthal_notch_phase);

  // Compute solid angle
  m_solid_angle = 2. * M_PI * ( cos(M_PI*m_zenithmin/180.) - cos(M_PI*m_zenithmax/180.) ); // (sr)

  // Integral flux
  m_integral_flux = m_solid_angle * m_fSpectralEnergy.Integral(m_emin, m_emax); // (m^-2 s^-1)

  // Integral intensity
  double integral_intensity=m_fSpectralEnergy.Integral(m_emin, m_emax); // (m^-2 s^-1 sr^-1)

  // Setting precision for inverse cumulative distribution functions
  m_cdf_steps=1000;
  m_cdf_energy_slices=100;
  double cdf_steps=static_cast<double> (m_cdf_steps);
  double cdf_energy_slices=static_cast<double> (m_cdf_energy_slices);

  std::cout << "Creating inverse CDF functions for Earth limb model. " 
	    << m_cdf_steps << " CDF steps and " << m_cdf_energy_slices << " slices in energy." 
	    << std::endl;

  // Create energy inverse cumulative distribution function
  double sum=0.;
  double log10_energy=log10(m_emin);
  double d_log10_energy=(log10(m_emax)-log10(m_emin))/cdf_steps;

  for(int ii=0;ii<m_cdf_steps;ii++){
    sum+=m_fSpectralEnergy.Eval(pow(10.,log10_energy))*pow(10.,log10_energy)*log(10.)*d_log10_energy;
    m_energy_inverse_cdf.SetPoint(m_energy_inverse_cdf.GetN(),sum,log10_energy);
    log10_energy+=d_log10_energy;
  }
  normalize_cdf(m_energy_inverse_cdf);

  // Create zenith and azimuth inverse cumulative distribution functions
  double zenith;
  double d_zenith=(m_zenithmax-m_zenithmin)/cdf_steps;
  double integral_zenith;
  double azimuth;
  double d_azimuth=(m_azimuthmax-m_azimuthmin)/cdf_steps;
  double integral_azimuth;

  log10_energy=log10(m_emin);
  d_log10_energy=(log10(m_emax)-log10(m_emin))/cdf_energy_slices;
  for(int ii=0;ii<m_cdf_energy_slices;ii++){

    // Set energy-dependent zenith parameters
    m_fZenith.SetParameter(2,m_fZenithSlope.Eval(log10_energy));
    m_fZenith.SetParameter(3,m_fZenithCurve.Eval(log10_energy));
    
    integral_zenith=m_fZenith.Integral(m_zenithmin,m_zenithmax);

    // Set energy-dependent azimuth parameters
    m_fAzimuth.SetParameter(2,m_fAzimuthalLogsine.Eval(log10_energy));
    m_fAzimuth.SetParameter(3,m_fAzimuthalNotchRatio.Eval(log10_energy));
    m_fAzimuth.SetParameter(4,m_fAzimuthalNotchWidth.Eval(log10_energy));

    integral_azimuth=m_fAzimuth.Integral(m_azimuthmin,m_azimuthmax);

    m_zenith_inverse_cdf.push_back(TGraph());
    m_azimuth_inverse_cdf.push_back(TGraph());

    // Create zenith inverse cumulative distribution function
    zenith=m_zenithmin;
    sum=0.;
    for(int jj=0;jj<m_cdf_steps;jj++){
      sum+=m_fZenith.Eval(zenith)*d_zenith;
      m_zenith_inverse_cdf[ii].SetPoint(m_zenith_inverse_cdf[ii].GetN(),sum/integral_zenith,zenith);
      zenith+=d_zenith;
    }
    normalize_cdf(m_zenith_inverse_cdf[ii]);
      
    // Create azimuth inverse cumulative distribution function
    azimuth=m_azimuthmin;
    sum=0.;
    for(int jj=0;jj<m_cdf_steps;jj++){
      sum+=m_fAzimuth.Eval(azimuth)*d_azimuth;
      m_azimuth_inverse_cdf[ii].SetPoint(m_azimuth_inverse_cdf[ii].GetN(),sum/integral_azimuth,azimuth);
      azimuth+=d_azimuth;
    }
    normalize_cdf(m_azimuth_inverse_cdf[ii]);

    // Output the azimuthal cumulative distribution function (debugging)
    //for(int jj=0;jj<m_cdf_steps;jj++){
    //  double x, y;
    //  m_azimuth_inverse_cdf[ii].GetPoint(jj, x, y);
    //  std::cerr << log10_energy << " " << x << " " << y << std::endl;
    //}

    log10_energy+=d_log10_energy;
  }

}

void EarthPhenomLimb::calculate(double &zenith, double &azimuth, double &energy) {
  double temp_zenith, temp_azimuth, temp_energy;

  // Calculate zenith angle (deg), azimuth (deg), and energy (MeV)
  double r_energy=CLHEP::RandFlat::shoot(),
         r_zenith=CLHEP::RandFlat::shoot(),
         r_azimuth=CLHEP::RandFlat::shoot();

  temp_energy=pow(10.,m_energy_inverse_cdf.Eval(r_energy));

  double cdf_energy_slices=static_cast<double> (m_cdf_energy_slices);
  int log10_energy_index=static_cast<int> ((((log10(temp_energy)-log10(m_emin))/(log10(m_emax)-log10(m_emin)))*cdf_energy_slices)+0.5);

  temp_zenith=m_zenith_inverse_cdf[log10_energy_index].Eval(r_zenith);
  temp_azimuth=m_azimuth_inverse_cdf[log10_energy_index].Eval(r_azimuth);
  
  // Invert photon directions (i.e. turn the Earth inside out for simulating back-entering events)
  if(m_invert_direction){
    zenith = 180.-temp_zenith;
    if(temp_azimuth<180.)
      azimuth = temp_azimuth+180.;
    else
      azimuth = temp_azimuth-180.;
  }else{
    zenith = temp_zenith;
    azimuth = temp_azimuth;
  }
  energy = temp_energy;

  return;
}

double EarthPhenomLimb::energy(double time) {
  (void)(time); // Earth limb emission is not variable in this version
  calculate(m_zenith,m_azimuth,m_energy);
  m_energy_called = true;
  return m_energy; // (MeV)
}

std::pair<double, double> EarthPhenomLimb::dir(double energy) {
  if(energy != m_energy || !m_energy_called){
    std::cerr << "ERROR in routine calling EarthPhenomLimb: need to call energy() before dir()." << std::endl;
    throw std::runtime_error("EarthPhenomLimb: ERROR in routine calling EarthPhenomLimb: need to call energy() before dir().");
  }
  m_energy_called = false;

  // Using zenith-local coordinates (cos(zenith), azimuth)
  // Azimuth is 0 in the north and 90 deg in the east when looking at the Earth
  return std::make_pair(cos(M_PI*m_zenith/180.), M_PI*m_azimuth/180.); // Convert from degrees to radians
}

std::pair<double, std::pair<double, double> > EarthPhenomLimb::photon(){
  calculate(m_zenith,m_azimuth,m_energy);

  // Using zenith-local coordinates (cos(zenith), azimuth)
  // Azimuth is 0 in the north and 90 deg in the east when looking at the Earth
  std::pair<double, double> direction = std::make_pair(cos(M_PI*m_zenith/180.), M_PI*m_azimuth/180.); // Convert from degrees to radians
  return std::make_pair(m_energy, direction); // Energy (MeV) 
}

ISpectrumFactory &EarthPhenomLimbFactory() {
   static SpectrumFactory<EarthPhenomLimb> myFactory;
   return myFactory;
}

void EarthPhenomLimb::normalize_cdf(TGraph &g){
  double a,b;
  g.GetPoint(g.GetN()-1,a,b);
  double cdf_normalization=a;
  for(int ii=0;ii<m_cdf_steps;ii++){
    g.GetPoint(ii,a,b);
    g.SetPoint(ii,a/cdf_normalization,b);
  }
}
