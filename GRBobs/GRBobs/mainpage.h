/*! \mainpage GRBobs Package: phenomenological model for GRB
   
This package is a redefinition of the GRB phenomenological model for Gamma-Ray Bursts

The simulator uses the same schema as the GRB physical model (GRB package), but, 
instead of exploring the physics for deriving light curves and the spectra, it uses phenomenological prescriptions.
The reference articles are listed below:

1) D. Band,et al. "BATSE observations of gamma-ray burst spectra. I Spectral diversity." Ap. J., 413:281-292, August 1993.

2) R. D. Preece, et al. "The BATSE Gamma-Ray Burst Spectral Catalog. I. High Time Resolution Spectroscopy of Bright Bursts Using High Energy Resolution Data." Ap. J. Supp., 126:19-36, January 2000.

3) E. E. Fenimore, et al. "Gamma-Ray Burst Peak Duration as a Function of Energy". Ap. J. Lett., 448:L101+, August 1995.

4) J. P. Norris, et al. "Attributes of Pulses in Long Bright Gamma-Ray Bursts." Ap. J., 459:393-+, March 1996.

The xml spectrum object has to be defined like:

\verbatim 
    <source name=" GRB1 ">
        <spectrum escale="MeV"> 
        <SpectrumClass name="GRBobsmanager" params="100,0.0e-5,20,-1.0,-2.25,100, 1"/>
        <celestial_dir ra="10." dec="22."/>
    </spectrum> </source>
\endverbatim
where the params are:
- Starting time of the burst (s)
- Fluence in \f$ erg/cm^{2}\f$. If it is 0 then it has sampled from the BATSE fluence distribution.
- The number of pulses of the bursts.
- The low energy spectral index (\f$\alpha\f$)
- The high energy spectral index (\f$\beta\f$)
- The minimim energy for extractint photons (in MeV).
- The flag for the GBM output. (1==YES and 0=NO). If GBM output is generated, the sumulator is slower and additional text file will be created 
(see, GRBobsSim::SaveGBMDefinition and GRBobsSim::GetGBMFlux)

<br>
The file  <a href="../../src/test/GRBParam.txt"> GRBParam.txt </a> holds the parameters for the ROOT test program.
<br>
It MUST be in /src/test directory of the GRBobs package. It contains: 
<br>
- The GRB number for initializing the Random number generator. Equal seeds are equal sequence of numbers
- The galacic l and b coordinates. NOTE if l or b are < -180 than their value is extracted randomly.
- The fluence in BATSE energy range (20 keV-1 MeV). NOTE if the value is zero than the fluence is set by Parameters::GetBATSEFluence()
- The Number of peaks to generate
- The characteristic time scale for extracting the interval between pulses.

<br>
\section test How-to use the test programs


The GRB simulator can be use with several test program:

- test_GRBROOT.exe
This test program makes use of the \e ROOT graphical environment to display some plots regarding the simulated GRB. It shows the evolution of the flux with respect to time, and plots the integrated spectrum and the light curves.
to execute it type "test_GRBROOT.exe".
Type "test_GRBROOT.exe -help" to have a list of the available options.
<br>

<hr>
\section requirements CMT requirements
\include cmt/requirements
<hr>
\section notes Release Notes
\include doc/release.notes
<hr>
 */


