/**
* @file EarthPhenomInner.cxx
* @brief A phenomenological model of the inner Earth based on LAT measurements
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

#include "EarthPhenom/EarthPhenomInner.h"
#include "TString.h"
#include "TGraph.h"

#include "CLHEP/Random/RandFlat.h"

//static SpectrumFactory<EarthPhenomInner> factory;
//const ISpectrumFactory& EarthPhenomInnerFactory = factory;

EarthPhenomInner::EarthPhenomInner(const std::string &paramString) {	
  std::vector<std::string> params;
  facilities::Util::stringTokenize(paramString, ", ", params);

  m_normalization = std::atof(params[0].c_str()); // Normalization = 1 is default
  m_emin = std::atof(params[1].c_str()); // MeV
  m_emax = std::atof(params[2].c_str()); // MeV
  if(std::atoi(params[3].c_str())!=0) // 0 = normal mode, 1 = invert direction
    m_invert_direction = true;
  else
    m_invert_direction = false;

  if(m_emax>2000.)
    std::cerr << "WARNING: maximum energy = " << m_emax 
	      << " MeV. Suggest maximum energy of 1000 MeV (may encounter numerical issues otherwise). Continue at your own risk..." << std::endl; 

  init(m_normalization, m_emin, m_emax);

  std::cerr << "EarthPhenomInner created. Normalization = "
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

double EarthPhenomInner::flux(double time) const { // Argument is the mission elapsed time (s)
    return m_integral_flux / solidAngle(); // (m^-2 s^-1 sr^-1)
}

double EarthPhenomInner::solidAngle() const {
  return m_solid_angle; // (sr)
}

void EarthPhenomInner::init(double normalization, double emin, double emax) {
  
  // Initialize pre-fitted model parameters

  // Spectral component
  m_spectral_prefactor = normalization * 1.e4 * 1.724e-02; // (MeV^-1 m^-2 s^-1 sr^-1) Converting from cm^-2 to m^-2
  m_spectral_index1 = -1.533;
  m_spectral_index2 = -4.546e+02;
  m_spectral_ebreak = 5.759e+02; // (MeV)
  m_spectral_beta = 9.476e+01;

  // Zenith angle component
  // Zenith = 0 deg -> away from Earth center; zenith = 180 deg -> towards Earth center
  m_zenithmin = 180.-50.; // (deg)
  m_zenithmax = 180.-0.; // (deg)
  
  //m_zenith_energy_slope_pol0 = -8.484e-02;
  //m_zenith_energy_slope_pol1 = 3.095e-01;
  //m_zenith_energy_slope_pol2 = -2.145e-01;
  //m_zenith_energy_slope_pol3 = 4.445e-02;
  m_zenith_energy_slope_pol0 = -4.031e-02; // Updated fit
  m_zenith_energy_slope_pol1 = 2.437e-01;
  m_zenith_energy_slope_pol2 = -1.829e-01;
  m_zenith_energy_slope_pol3 = 3.950e-02;

  // Azimuth angle component
  // Looking towards the Earth
  // Simulation coordinate system: Azimuth = 0 deg -> east; azimuth = 90 deg -> north (right-handed)
  // LAT analysis coordinate system: Azimuth = 0 deg -> north; azimuth = 90 deg -> east (left-handed)
  m_azimuthmin = 0.; // (deg)
  m_azimuthmax = 360.; // (deg)
  m_azimuthal_logsine_phase = 90.-176.5; // (deg) Converting from LAT analysis to simulation coordinate system

  //m_azimuthal_energy_logsine_pol0 = 2.724e-01;
  //m_azimuthal_energy_logsine_pol1 = -2.373e-01;
  //m_azimuthal_energy_logsine_pol2 = 8.176e-02;
  //m_azimuthal_energy_logsine_pol3 = 3.357e-02;
  //m_azimuthal_energy_logsine_pol4 = -1.244e-02;
  m_azimuthal_energy_logsine_pol0 = 2.762e-01; // Updated fit
  m_azimuthal_energy_logsine_pol1 = -2.427e-01;
  m_azimuthal_energy_logsine_pol2 = 8.318e-02;
  m_azimuthal_energy_logsine_pol3 = 3.425e-02;
  m_azimuthal_energy_logsine_pol4 = -1.265e-02;

  // Initialize model formulae

  // Variable is energy (MeV, linear-space)
  TString spectralEnergyFormula = "(x<(10^2.6))*[0]*x^[1]*(1+(x/[3])^(([1]-[2])/[4]))^(-1.*[4])";
  m_fSpectralEnergy = TF1("fSpectralEnergy", spectralEnergyFormula.Data(), m_emin, m_emax);
  m_fSpectralEnergy.SetParameter(0,m_spectral_prefactor);
  m_fSpectralEnergy.SetParameter(1,m_spectral_index1);
  m_fSpectralEnergy.SetParameter(2,m_spectral_index2);
  m_fSpectralEnergy.SetParameter(3,m_spectral_ebreak);
  m_fSpectralEnergy.SetParameter(4,m_spectral_beta);

  // Variable is energy (MeV, log10-space)
  TString spectralFormula = "(x<2.6)*[0]*(10^x)^[1]*(1+((10^x)/[3])^(([1]-[2])/[4]))^(-1.*[4])"; // Weighted for probability
  m_fSpectral = TF1("fSpectral", spectralFormula.Data(), log10(m_emin), log10(m_emax));
  m_fSpectral.SetParameter(0,m_spectral_prefactor);
  m_fSpectral.SetParameter(1,m_spectral_index1);
  m_fSpectral.SetParameter(2,m_spectral_index2);
  m_fSpectral.SetParameter(3,m_spectral_ebreak);
  m_fSpectral.SetParameter(4,m_spectral_beta);

  // Variable is energy (MeV, log10-space)
  // Energy-dependent variation in inner Earth logarithmic slope
  TString zenithSlopeFormula = "[0]+[1]*x+[2]*x^2+[3]*x^3";
  m_fZenithSlope = TF1("fZenithSlope", zenithSlopeFormula.Data(), log10(m_emin), log10(m_emax));
  m_fZenithSlope.SetParameter(0,m_zenith_energy_slope_pol0);
  m_fZenithSlope.SetParameter(1,m_zenith_energy_slope_pol1);
  m_fZenithSlope.SetParameter(2,m_zenith_energy_slope_pol2);
  m_fZenithSlope.SetParameter(3,m_zenith_energy_slope_pol3);
  
  // Variable is zenith angle (deg)
  // Parameter 0 is energy-dependent
  // Notice that zenith angle is used instead of nadir angle (sign flip)
  TString zenithFormula = "exp([0]*(180.-x))";
  m_fZenith = TF1("fZenith", zenithFormula.Data(), m_zenithmin, m_zenithmax);
  //m_fZenith.SetParameter(0,m_zenith_slope);

  // Variable is energy (MeV, log10-space)
  // Energy-dependent variation of logsine component amplitude
  TString azimuthalLogsineFormula = "max([0]+[1]*x+[2]*x^2+[3]*x^3+[4]*x^4,0)";
  m_fAzimuthalLogsine = TF1("fAzimuthalLogsine", azimuthalLogsineFormula.Data(), log10(m_emin), log10(m_emax));
  m_fAzimuthalLogsine.SetParameter(0,m_azimuthal_energy_logsine_pol0);
  m_fAzimuthalLogsine.SetParameter(1,m_azimuthal_energy_logsine_pol1);
  m_fAzimuthalLogsine.SetParameter(2,m_azimuthal_energy_logsine_pol2);
  m_fAzimuthalLogsine.SetParameter(3,m_azimuthal_energy_logsine_pol3);
  m_fAzimuthalLogsine.SetParameter(4,m_azimuthal_energy_logsine_pol4);

  // Variable is azimuth angle (deg)
  // Parameter 1 is energy-dependent
  // The "-" sign in front of azimuth angle converts from LAT analysis to simulation coordinate system
  TString azimuthalFormula = "exp([1]*sin(TMath::Pi()*(-x-[0])/180))";
  m_fAzimuth = TF1("fAzimuth", azimuthalFormula.Data(), m_azimuthmin, m_azimuthmax);
  m_fAzimuth.SetParameter(0,m_azimuthal_logsine_phase);

  // Compute solid angle
  m_solid_angle = 2. * M_PI * ( cos(M_PI*m_zenithmin/180.) - cos(M_PI*m_zenithmax/180.) ); // (sr)

  // Integral flux
  m_integral_flux = m_solid_angle * m_fSpectralEnergy.Integral(m_emin, m_emax); // (m^-2 s^-1)

  // Setting precision for inverse cumulative distribution functions
  m_cdf_steps=1000;
  m_cdf_energy_slices=100;
  double cdf_steps=static_cast<double> (m_cdf_steps);
  double cdf_energy_slices=static_cast<double> (m_cdf_energy_slices);

  std::cout << "Creating inverse CDF functions for inner Earth model. " 
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

    // Set energy-dependent zenith parameter
    m_fZenith.SetParameter(0,m_fZenithSlope.Eval(log10_energy));

    // Set energy-dependent azimuth parameter
    m_fAzimuth.SetParameter(1,m_fAzimuthalLogsine.Eval(log10_energy));
    
    //integral_azimuth=m_fAzimuth.Integral(m_azimuthmin,m_azimuthmax);

    m_zenith_inverse_cdf.push_back(TGraph());
    m_azimuth_inverse_cdf.push_back(TGraph());

    // Create zenith inverse cumulative distribution function
    zenith=m_zenithmin;
    sum=0.;
    for(int jj=0;jj<m_cdf_steps;jj++){
      sum+=m_fZenith.Eval(zenith)*(cos(M_PI*zenith/180.)-cos(M_PI*(zenith+d_zenith)/180.));
      m_zenith_inverse_cdf[ii].SetPoint(m_zenith_inverse_cdf[ii].GetN(),sum,zenith);
      zenith+=d_zenith;
    }
    normalize_cdf(m_zenith_inverse_cdf[ii]);

    // Create azimuth inverse cumulative distribution function
    azimuth=m_azimuthmin;
    sum=0.;
    for(int jj=0;jj<m_cdf_steps;jj++){
      sum+=m_fAzimuth.Eval(azimuth)*d_azimuth;
      //m_azimuth_inverse_cdf[ii].SetPoint(m_azimuth_inverse_cdf[ii].GetN(),sum/integral_azimuth,azimuth);
      m_azimuth_inverse_cdf[ii].SetPoint(m_azimuth_inverse_cdf[ii].GetN(),sum,azimuth);
      azimuth+=d_azimuth;
    }
    normalize_cdf(m_azimuth_inverse_cdf[ii]);

    log10_energy+=d_log10_energy;
  }

}

void EarthPhenomInner::calculate(double &zenith, double &azimuth, double &energy) {
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

double EarthPhenomInner::energy(double time) {
  (void)(time); // Earth limb emission is not variable in this version
  calculate(m_zenith,m_azimuth,m_energy);
  m_energy_called = true;
  return m_energy; // (MeV)
}

std::pair<double, double> EarthPhenomInner::dir(double energy) {
  if(energy != m_energy || !m_energy_called){
    std::cerr << "ERROR in routine calling EarthPhenomInner: need to call energy() before dir()." << std::endl;
    throw std::runtime_error("EarthPhenomInner: ERROR in routine calling EarthPhenomInner: need to call energy() before dir().");
  }
  m_energy_called = false;

  // Using zenith-local coordinates (cos(zenith), azimuth)
  // Azimuth is 0 in the east and 90 deg in the north when looking at the Earth
  return std::make_pair(cos(M_PI*m_zenith/180.), M_PI*m_azimuth/180.); // Convert from degrees to radians
}

std::pair<double, std::pair<double, double> > EarthPhenomInner::photon(){
  calculate(m_zenith,m_azimuth,m_energy);

  // Using zenith-local coordinates (cos(zenith), azimuth)
  // Azimuth is 0 in the east and 90 deg in the north when looking at the Earth
  std::pair<double, double> direction = std::make_pair(cos(M_PI*m_zenith/180.), M_PI*m_azimuth/180.); // Convert from degrees to radians
  return std::make_pair(m_energy, direction); // Energy (MeV) 
}

ISpectrumFactory & EarthPhenomInnerFactory() {
   static SpectrumFactory<EarthPhenomInner> myFactory;
   return myFactory;
}

void EarthPhenomInner::normalize_cdf(TGraph &g){
  double a,b;
  g.GetPoint(g.GetN()-1,a,b);
  double cdf_normalization=a;
  for(int ii=0;ii<m_cdf_steps;ii++){
    g.GetPoint(ii,a,b);
    g.SetPoint(ii,a/cdf_normalization,b);
  }
}
