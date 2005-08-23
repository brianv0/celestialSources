// Mainpage for doxygen

/** @mainpage package genericSources

 @author James Chiang

 @section intro Introduction

 The sources in this package implement flux::ISpectrum sources to
 provide photons from astrophysical sources for Gleam and
 observationSim.

 <hr>
 @section notes Release Notes
  release.notes

 <hr>
 @section requirements requirements
 @verbinclude requirements
*/

/**
 @page sources genericSources

 @section available_sources Source Definitions

 In addition to the sources that are already available from the flux
 package (e.g., constant point sources with power-law or broken
 power-law spectra, isotropic diffuse in instrument coordinates,
 MapSpectrum's Galactic diffuse), there are many sources that are
 available to gtobssim and Gleam through the celestialSources
 container package. Here we document the sources available in the
 genericSources subpackage. These sources derive from the Spectrum
 class of the flux package, and accordingly, their xml entries all
 have the same form. The params string circumvents the DTD and thereby
 allows essentially free-format data to be passed to the Spectrum
 class.

 For each source, we describe the entries in the <tt>param</tt>
 string.  Optional parameters have their default values given in
 parentheses.  

 - GaussianSource
   - <b>flux</b> Total flux in units of \f$\mbox{m}^{-2}\mbox{s}^{-1}\f$.
   - <b>gamma</b> Photon spectral index such that 
      \f$dN/dE \propto E^{-\Gamma}\f$.
   - <b>RA</b> Source centroid J2000 right ascension in degrees.
   - <b>Dec</b> Source centroid J2000 declination in degrees.
   - <b>major axis (1)</b> Semi-major axis of the 68\% CL contour in degrees.
   - <b>minor axis (1)</b> Semi-minor axis of the 68\% CL contour in degrees.
   - <b>position angle (0)</b> Position angle of the major axis measured 
     wrt North in degrees.
   - <b>Emin (30)</b> Minimum photon energy in MeV.
   - <b>Emax (1e5)</b> Maximum photon energy in MeV.
@verbatim
   <source name="gaussian_source">
      <spectrum escale="MeV">
         <SpectrumClass name="GaussianSource"
                        params="0.1, 2.1, 45., 30., 3., 0.5, 45, 30., 2e5"/>
         <use_spectrum frame="galaxy"/>
      </spectrum>
   </source>
@endverbatim

 - Isotropic
   - <b>flux</b> in units of \f$\mbox{m}^{-2}\mbox{s}^{-1}\mbox{sr}^{-1}\f$.
   - <b>gamma</b> Photon spectral index such that 
      \f$dN/dE \propto E^{-\Gamma}\f$.
   - <b>Emin (30)</b> Minimum photon energy in MeV.
   - <b>Emax (1e5)</b> Maximum photon energy in MeV.
@verbatim
   <source name="Extragalactic_diffuse">
      <spectrum escale="MeV">
         <SpectrumClass name="Isotropic"
                        params="10.7, 2.1, 20., 2e5"/>
         <use_spectrum frame="galaxy"/>
      </spectrum>
   </source>
@endverbatim

 - MapCube
   - <b>flux</b> Total flux from the map, integrated over solid angle, in 
     units of \f$\mbox{m}^{-2}\mbox{s}^{-1}\f$.
   - <b>FITS file</b> A plate-carree FITS image in Galactic or J2000 
     coordinates.
@verbatim
   <source name="map_cube_source">
      <spectrum escale="MeV">
         <SpectrumClass name="MapCube"
          params="1., $(GENERICSOURCESROOT)/data/test_image.fits"/>
         <use_spectrum frame="galaxy"/>
      </spectrum>
   </source>
@endverbatim

 - MapSource
   - <b>flux</b> Total flux from the map, integrated over solid angle, in 
     units of \f$\mbox{m}^{-2}\mbox{s}^{-1}\f$.
   - <b>gamma</b> Photon spectral index such that 
      \f$dN/dE \propto E^{-\Gamma}\f$.
   - <b>FITS file</b> A plate-carree FITS image in Galactic or J2000 
     coordinates.
   - <b>Emin (30)</b> Minimum photon energy in MeV.
   - <b>Emax (1e5)</b> Maximum photon energy in MeV.
@verbatim
<!-- MapSource version of the Galactic Diffuse model -->
   <source name="Galactic_diffuse">
      <spectrum escale="MeV">
         <SpectrumClass name="MapSource"
          params="17.,2.1,$(FLUXROOT)/sources/gas_gal.fits,30.,2e5"/>
         <use_spectrum frame="galaxy"/>
      </spectrum>
   </source>
@endverbatim

 - PeriodicSource A point source with sinusoidal light curve.
   - <b>flux</b> Average flux in units of \f$\mbox{m}^{-2}\mbox{s}^{-1}\f$.
   - <b>gamma</b> Photon spectral index such that 
      \f$dN/dE \propto E^{-\Gamma}\f$.
   - <b>period</b> Source period in seconds.
   - <b>amplitude (0.5)</b> Amplitude of the sinusoidal modulation.
   - <b>phi0 (0)</b> Phase offset specified on the unit interval.
   - <b>Emin (30)</b> Minimum photon energy in MeV.
   - <b>Emax (1e5)</b> Maximum photon energy in MeV.
@verbatim
   <source name="periodic_source">
      <spectrum escale="MeV">
         <SpectrumClass name="PeriodicSource"
                        params="0.1, 2.1, 1e3, 1, 0.75, 30., 2e5"/>
         <galactic_dir l="0" b="0"/>
      </spectrum>
   </source>
@endverbatim

 - Pulsar A pulsar source whose light curve is given by an ascii template
   file.
   - <b>flux</b> Average flux in units of \f$\mbox{m}^{-2}\mbox{s}^{-1}\f$.
   - <b>gamma</b> Photon spectral index such that 
      \f$dN/dE \propto E^{-\Gamma}\f$.
   - <b>period</b> Pulsar period in seconds.
   - <b>pdot</b> Time derivative of the pulsar period in 
     \f$\mbox{s}\,\mbox{s}^{-1}\f$.
   - <b>t0</b> Reference epoch in MET seconds.
   - <b>template file</b> Filename of the ascii light curve template.
     The file should consist of two columns, phase and intensity.  The
     phase intervals must be uniformly spaced.  The phase scale and
     absolute intensities are arbitrary and rescaled using the flux
     and period parameter values.
   - <b>phi0 (0)</b> Phase offset in the unit interval.
   - <b>Emin (30)</b> Minimum photon energy in MeV.
   - <b>Emax (1e5)</b> Maximum photon energy in MeV.
@verbatim
   <source name="Crab_Pulsar">
      <spectrum escale="MeV">
         <SpectrumClass name="Pulsar"
         params="1e-3,2.,0.033,0,0,$(OBSERVATIONSIMROOT)/data/CrabTemplate.dat"/>
         <celestial_dir ra="83.57" dec="22.01"/>
      </spectrum>
   </source>
@endverbatim

 - SimpleTransient A point source with a single active interval during
   which it has a constant flux and power-law spectrum.
   - <b>flux</b> Flux while in the active state in units of 
     \f$\mbox{m}^{-2}\mbox{s}^{-1}\f$.
   - <b>gamma</b> Photon spectral index such that 
      \f$dN/dE \propto E^{-\Gamma}\f$.
   - <b>tstart</b> Start time of the active state in MET seconds.
   - <b>tstop</b> Stop time of the active state in MET seconds.
   - <b>Emin (30)</b> Minimum photon energy in MeV.
   - <b>Emax (1e5)</b> Maximum photon energy in MeV.
@verbatim
   <source name="simple_transient">
      <spectrum escale="MeV">
         <SpectrumClass name="SimpleTransient"
          params="10., 2., 1e3, 1.1e3, 30., 2e5"/>
         <celestial_dir ra="83." dec="22."/>
      </spectrum>
   </source>
@endverbatim

 - SpectralTransient
   - <b>flux</b> Mean flux during the active state in units of 
     \f$\mbox{m}^{-2}\mbox{s}^{-1}\f$.
   - <b>tstart</b> Start time of the active state in MET seconds.
   - <b>tstop</b> Stop time of the active state in MET seconds.
   - <b>template file</b> Filename of the light curve template.  
      May be ascii or FITS.
   - <b>Emin (20)</b> Minimum photon energy in MeV.
   - <b>Emax (2e5)</b> Maximum photon energy in MeV.
   - <b>lc (0)</b> light curve number, if FITS file.
   - <b>z (0)</b> Redshift used for EBL attenuation calculation.
   - <b>useLogParabola (0)</b> Flag to use log-parabolic form for the
     spectrum rather than a broken power-law
@verbatim
   <source name="spectral_transient">
      <spectrum escale="MeV">
          <SpectrumClass name="SpectralTransient"
          params="flux=1e-1, tstart=0., tstop=1e4, templateFile=$(GENERICSOURCESROOT)/data/testTemplate.dat, emin=20, emax=2e5, lc=0, z=0, useLogParabola=0"/>
          <celestial_dir ra="193.4" dec="-5.82"/>
      </spectrum>
   </source>
@endverbatim

 - TransientTemplate A point source transient whose active state light curve
   shape is given by an ascii template file.
   - <b>flux</b> Mean flux during the active state in units of 
     \f$\mbox{m}^{-2}\mbox{s}^{-1}\f$.
   - <b>gamma</b> Photon spectral index such that 
      \f$dN/dE \propto E^{-\Gamma}\f$.
   - <b>tstart</b> Start time of the active state in MET seconds.
   - <b>tstop</b> Stop time of the active state in MET seconds.
   - <b>template file</b> Filename of the ascii light curve template.
     The file should consist of two columns, time and intensity.  The
     time intervals must be uniformly spaced.  The time scale and
     absolute intensities are arbitrary and rescaled using the flux,
     tstop, and tstep parameter values.
   - <b>Emin (30)</b> Minimum photon energy in MeV.
   - <b>Emax (1e5)</b> Maximum photon energy in MeV.
@verbatim
   <source name="transient_template">
      <spectrum escale="MeV">
         <SpectrumClass name="TransientTemplate"
          params="100.,2,1e3,1.1e3,$(OBSERVATIONSIMROOT)/data/CrabTemplate.dat"/>
         <celestial_dir ra="80" dec="20"/>
      </spectrum>
   </source>
@endverbatim

*/
