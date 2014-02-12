#ifndef mySpectrum_MapSource_h
#define mySpectrum_MapSource_h

#include "flux/Spectrum.h"
#include <vector>
#include <stdlib.h>
#include "TF1.h"
#include "TGraph.h"

class EarthPhenomLimb : public Spectrum {
public:

  /// This constructor is required and used in FluxSource for
  /// "SpectrumClass" sources.
  EarthPhenomLimb(const std::string &params);

  virtual ~EarthPhenomLimb(){}
	
  /// @return Particle type, "gamma".
  virtual const char * particleName() const {return "gamma";}
  
  /// @return average differential flux integrated over energy (m^-2 s^-1 sr^-1) time??? -> check this
  /// @param time Simulation time in seconds.
  virtual double flux(double time) const;

  /// @return "Effective" solid angle (sr).
  virtual double solidAngle() const;

  /// @return Title describing the spectrum.
  virtual std::string title() const {return "EarthPhenomLimb";}

  /// @return Interval to the next event (seconds)
  virtual double interval(double time){return -1;}; // flag for Poisson 

  /// @return Photon energy (MeV).
  virtual double energy(double time);

  /// @return Photon direction in (zenith angle (degrees), azimuth (degrees,East=0,North=90deg)).
  virtual std::pair<double, double> dir(double energy);

  /// @return energy and photon direction (zenith angle (degrees), azimuth (degrees,East=0,North=90deg)).
  std::pair<double, std::pair<double, double> > photon();
	 
  /// Initialize energy range (MeV)
  void init(double normalization, double emin, double emax);

private:
  /// Normalization, 1 should be default
  double m_normalization;

  /// User-specified energy range
  double m_emin, m_emax; // MeV, MeV

  // Invert photon directions, False = normal mode, True = invert photon directions
  bool m_invert_direction;

  /// For testing
  //bool m_simple_zenith;
  //bool m_simple_azimuth;

  /// Solid angle and integral flux of Earth limb, calculated in init()
  double m_solid_angle; // (sr)
  double m_integral_flux; // (m^-2 s^-1 sr^-1)

  /// Model formulae 
  TF1 m_fSpectral;
  TF1 m_fSpectralEnergy;
  TF1 m_fZenithSlope;
  TF1 m_fZenithCurve;
  TF1 m_fZenithBreak;
  TF1 m_fZenith;
  TF1 m_fAzimuthalLogsine;
  TF1 m_fAzimuthalNotchRatio;
  TF1 m_fAzimuthalNotchWidth;
  TF1 m_fAzimuth;

  /// Pre-fitted model parameters
  double m_spectral_prefactor;
  double m_spectral_index1;
  double m_spectral_index2;
  double m_spectral_ebreak;
  double m_spectral_beta;

  double m_zenithmin;
  double m_zenithmax;
  double m_zenith_peak;
  double m_zenith_width;

  double m_zenith_nadir_peak;
  double m_zenith_nadir_width;
  
  double m_zenith_theta_energy_slope_break_energy;
  double m_zenith_theta_energy_slope_le_prefactor;
  double m_zenith_theta_energy_slope_le_index;
  double m_zenith_theta_energy_slope_he_plateau_energy;
  double m_zenith_theta_energy_slope_he_plateau_halfwidth;
  
  double m_zenith_theta_energy_curve_he_plateau;
  double m_zenith_theta_energy_curve_prefactor;
  double m_zenith_theta_energy_curve_tanh_center;
  double m_zenith_theta_energy_curve_tanh_width;
  double m_zenith_theta_energy_curve_index;
  double m_zenith_theta_energy_curve_pol1;
  double m_zenith_theta_energy_curve_pol2;
  double m_zenith_theta_energy_curve_pol3;
  
  double m_zenith_theta_energy_break_prefactor;
  double m_zenith_theta_energy_break_tanh_center;
  double m_zenith_theta_energy_break_tanh_width;
  
  double m_azimuthmin;
  double m_azimuthmax;
  double m_azimuthal_logsine_phase;
  double m_azimuthal_notch_phase;

  double m_azimuthal_energy_logsine_prefactor;
  double m_azimuthal_energy_logsine_index;
  double m_azimuthal_energy_logsine_ecutoff;

  double m_azimuthal_energy_notch_ratio_prefactor;
  double m_azimuthal_energy_notch_ratio_index;
  double m_azimuthal_energy_notch_ratio_ecutoff;
  double m_azimuthal_energy_notch_width_prefactor;
  double m_azimuthal_energy_notch_width_index;
  double m_azimuthal_energy_notch_width_ecutoff;

  /// Inverse cumulative distribution functions

  int m_cdf_steps;
  int m_cdf_energy_slices;

  TGraph m_energy_inverse_cdf;
  std::vector<TGraph> m_zenith_inverse_cdf;
  std::vector<TGraph> m_azimuth_inverse_cdf;

  /// Internal helper function
  void calculate(double &zenith, double &azimuth, double &energy); // return zenith angle (deg), azimuth (deg), and energy (MeV)
  double m_zenith, m_azimuth, m_energy;
  double m_energy_called;

  /// Normalize CDF
  void normalize_cdf(TGraph &g);
};

#endif
