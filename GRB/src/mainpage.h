 /*! \mainpage GRB Package: simulating a transient source
   
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
 -# The phenomenological simulator of GRB signal. 
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
oriented magnetic field. The electrons can loose their energy via synchrotron
 emission. The characteristic synchrotron spectrum is boosted (and beamed) 
thanks to the Lorentz factor of the emitting material. 
The higher energy part of a GRB spectrum can be obtained keeping into account 
the possibility of Compton scattering of the synchrotron photons against the 
electron accelerated by the shock (Inverse Compton Scattering).

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
Depending on how GRBengine works, the sequence of shock that give up the GRB
can be formed in different ways. 
The File <a href="../../src/test/GRBParam.txt"> GRBParam.txt </a> contains information about:
- The GRB number for initializing the Random number generator. Equal seeds are equal sequence of numbers
- The galacic l and b coordinates. NOTE if l or b are < -180 than their value is extracted randomly.
- The fluence in BATSE energy range (20 keV-1 MeV). NOTE if the value is zero than the fluence is set by Parameters::GetBATSEFluence()
- The Number of shell
- The variability of the central engine which determines the characteristic temporal scale.
- The ratio between the maximum lorentz factor and the minimum lorentz factor of the shells
- The energy at which \f$\nu f(\nu)\f$ peaks (approximately)
- The ratio between the synchrotron and Inverse Compton energies (0-> Only synchrotron, 1 equal energy in synchrotron and in Ic);

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


