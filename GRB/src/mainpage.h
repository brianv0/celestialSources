 /*! \mainpage The GRB Physical model: simulating a transient source
   
The GRB package has been developed to provide the Glast simulator software 
with a full fledged simulator of a physical transient source.
As such it has been designed to correctly interface FluxSvc, 
by providing spectrum classes derived from ISpectrum.
This package is also developed to study the physics of bursts and the 
capability of Glast in the observation of rapid 
transient signals. It is interfaced to ROOT for possible visualization.
 
 This package contains four parts related to Gamma-Ray Bursts:
 -# The full physical simulator of GRB, based on the fireball model in 
 the internal shocks configuration.
 -# The original phenomenological simulator of GRB signal. 
 -# The algorithm describing the LAT alert.
 
 \section physics Very short introduction to the physical model
<br>
We adopted a physical model based on the \e fireball model of Gamma Ray 
Burst, because it is able to well reproduce the fast time variability 
observed in the GRB signal: a series of shells is injected in the circum 
burst medium with different Lorentz factors. 
When a faster shell reaches a slower one a shock occurs, and an accelerated 
electron distribution is obtained due to the shock acceleration process.
Some of the energy dissipated during the shock is converted into a randomly 
oriented magnetic field. The electrons can lose their energy via synchrotron
 emission. The characteristic synchrotron spectrum is boosted (and beamed) 
thanks to the Lorentz factor of the emitting material.
The doppler beaming determines, in this model, the pulse shape.
The higher energy part of a GRB spectrum is obtained by Inverse Compton scattering 
the synchrotron photons against the electron accelerated by the shock (Self Synchrotron Compton).

A short list of articles regarding the physics of GRB, 
observations and theories can be found <a href="http://www.pi.infn.it/~omodei/biblio.html"> here</a>. 

<br>

\subsection parameters Parameters for the GRB physical model
<br>
The file  <a href="../../src/test/GRBParam.txt"> GRBParam.txt </a> holds information used for the description of 
the GRB. 
It can be changed without recompiling the entire package.
<br>
It MUST be in /src/test directory of the GRB package. 
<br>
It contains the following observable quantities (typical values between parenthesis):
- Fluence between 20 keV and 1 MeV, in erg/cm^2 (BATSE values)
- Burst type: S=short bursts, L = long burst, R -> 25% S, 75% L
- Cut-off energy \f$ E_{co} \f$, in GeV. (3, 10 GeV) 
- Peak Energy \f$ E_{peak} \f$ of the synchrotron spectrum, in keV 
(it is sampled from a log-normal distribution with mean 235 keV and sigma 1.75. Moreover, \f$ E_{peak,SHORT}=E_{peak,LONG}/2 \f$)
- IC/Syn is the ratio between the Inverse Compton peak and the Synchrotron peak of the \f$ e^2 N(e)\f$ spectrum (0=pure synchrotron. < 10 typically)
- GBM flag, if 1 the program generates a series of files to be used with the GBM software
<br>
Notice: if the value of some parameter is set to 0, the Parameters class extract randomly the values from the appropriate distributions.
<br>
The variability scale \f$ t_v \f$ is set by the observations, requiring that for short bursts it is approximately equal to the burst duration 
(\f$ t_v\approx T_{90}\f$), whicle for long bursts it is sampled from the observed distribution of the FWHM: (\f$ t_v=10^{Gaus(0.0,0.5)}\f$).
The variability time scale (which is directly related to the pulse width) is contsrained to a minimum: 
The shortest variability time scale is a fraction (1/50) of the duration of the bursts: 
long bursts will likely have long pulses instead of many short pulses.
    
    This quantities are converted into the parameters of the model by a series of relations:
    
    - \f$ \Gamma=40.5 * E_{co} \f$
    - \f$ \frac{\Delta_0}{10^7}  =  13.6~(3~\alpha_B)(\frac{E_{tot}}{10^{52}})(\frac{3 \alpha_e \xi}{E_{co}})(\frac{E_{peak}}{100} t_v)^2 \f$
    - \f$ r_0  =  2~c~t_v \f$
    
    The total energy is set (and fixed) to  \f$ 10^{52} \f$ erg, and \f$ \Gamma_M/\Gamma_m=2 \f$, so that: 
    \f$ \Gamma_m = \frac{2~\Gamma}{\Gamma_M/\Gamma_m+1} \f$, and \f$ \Gamma_M= \frac{\Gamma_M}{\Gamma_m}\Gamma_m \f$ 
    This is a model dependent way to correlate observables and model parameters.

The parameters of the fireball are then computed in Parameters::ComputeParametersFromFile(std::string paramFile, int NGRB).

\section phenomenological Very short introduction to the phenomenological model

\section alert Very short introduction to the alert algorithm

<br>
\section test How-to use the test programs

The GRB simulator can be use with several test program:

- test_GRB.exe
This executable has been created for testing the code and it is a copy of the testMgr test program.
To launch it type "test_GRB.exe". The test program will load all the sources of the GRB_user_library.xml file.

- test_GRBROOT.exe
This test program makes use of the \e ROOT graphical environment to display some plots regarding the simulated GRB. It shows the evolution of the flux with respect to time, and plots the integrated spectrum and the light curves.
to execute it type "test_GRBROOT.exe".
Type "test_GRBROOT.exe -help" to have a list of the available options.
See <a href="../../src/test/other/GRBROOTtest.cxx"> GRBROOTtest.cxx</a> 
the test code.
<br>

<hr>
\section requirements CMT requirements
\include cmt/requirements
<hr>
\section notes Release Notes
\include doc/release.notes
<hr>
 */


