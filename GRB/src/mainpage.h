 /*! \mainpage GRB Package: simulating a transient source
   
The GRB package has been developed to provide the GLAST simulator software 
with a full fledged simulator of a physical transient source.
As such it has been designed to correctly interface FluxSvc, 
by providing spectrum classes derived from ISpectrum.
This package is also developed to study the physics of bursts and the 
capability of Glast in the observation of rapid 
transient signals. It is interfaced to ROOT for possible visualization.
 
 This package contains three parts related to Gamma-Ray Bursts:
 -# The full physical simulator of GRB, based on the fireball model in 
the internal shocks configuration.
 -# The phenomenological simulator of GRB signal. 
 -# The algorithm describing the LAT alert.

 
 \section physics Very short introduction to the physical model
<br>
We adopted a physical model based on the \em fireball model of Gamma Ray Burst, because it is able to well reproduce the fast time variability observed in the GRB signal: a series of shells is injected in the circum burst medium with different Lorentz factor. When a faster shell reach a slower one a shock occurs, and an accelerated electron distribution is obtained due to the shock acceleration process. Some of the energy dissipated during the shock is converted into a randomly oriented magnetic field. The electrons can loose their energy via synchrotron emission. The characteristic synchrotron spectrum is boosted (and beamed) thanks to the Lorentz factor of the emitting material. The higher energy part of a GRB spectrum can be obtained keeping into account the possibility of Compton scattering of the synchrotron photons against the electron accelerated by the shock (Inverse Compton Scattering).
<br>

\subsection parameters Parameters for the GRB physical model
<br>
The file /src/test/GRBParam.txt holds information used for the description of 
the GRB. 
It can be changed without recompiling the entire package.
<br>
It MUST be in /src/test directory of the GRB package. 
<br>
Depending on how GRBengine works, the sequence of shock that give up the GRB
can be formed in different ways. 
The File GRBParam.txt contains information about:
- The Number of shocks
- The total energy available to the fireball
- The Redshift of the source
- The minimum energy of the extracted photons. Only photons with a greater 
energy are extracted from the flux and stent to FluxSvc
- The shell type, jet or iso. See GRBShell for explanations.
- The radius of the jet (the jet is thought to be a cylinder)
- The angle between the jet direction and the line of sight
- The type of engine (see GRBengine for further explanations)
- The observed duration of the burst
- The observed rise time of the spikes
- The observed decay time of the spikes
- The energy at which \f$\nu f(\nu)\f$ peaks
- The thickness of the shells
- The minimum Lorentz factor to assign to a shell
- The maximum Lorentz factor to assign to a shell

A tipical "GRBParam.txt" configuration file look like this:
\verbatim
1	        // Nshock
1.0e+55		// Etot
1.0		// Redshift
1.0e+7          // Minimum Energy Extracted
Shell  type = 1 // 1= jet, 0 = spherical 
## If Shell == 1 : Jet Shells
1.0e+18	  	// Radius of the jet [cm]
0.0		// Angle of the jet  [degree]
## If Shell == 0 : Isotropic Shells
1.0e+18	  	// Radius of the Shell [cm]	
Engine type = 0 // 0=OP, 1= PP1S, 2= PP2S
## If  Engine type == 0 (Observed parameters are):
.2		// Duration of the Burst.
.1 		// Rise  time of a peak
.1              // decay time of a peak
0.3             // Peak  energy (MeV)
## If  Engine type == 1 or 2 (Physical parameters are):
1.0e+13		// Thickness of the shells[cm]
1000.		// Minimum Lorenz Factor
1000.		// Maximum Lorenz Factor
\endverbatim


\section phenomenological Very short introduction to the phenomenological model

\section alert Very short introduction to the alert algorithm

<br>
\section test How-to use the test programs


The GRB simulator can be use with several test program:

- test_GRB.exe
This executable tests the GRB algorithm. It initializes the GRB simulation,
and extracts photons according to the computed spectrum. 
To launch it type "test_GRB.exe". All the options are contained in the 
joboptions.txt file.

- test_GRBROOT.exe
This test program makes use of the \em ROOT graphical environment to display some plots regarding the simulated GRB. It shows the evolution of the flux with respect to time, and plots the integrated spectrum and the light curves.
to execute it type "test_GRBROOT.exe".
Type "test_GRBROOT.exe -help" to have a list of the available execution 
options.
<br>

\section Howto How-to use the GRB spectrum in Gleam
<br>
GRB is an independent package, it is an external service and it has to be 
declared in the joboption file of FluxSvc. For example, here are the 
two lines that add GRBSvc in the external services available and add the 
GRB shared library. 
They could be added, for example, to the "defaultOptions.txt" of FluxSvc
<br>
\verbatim 
//================================================== 
ApplicationMgr.ExtSvc += { "GRBSvc" };
ApplicationMgr.DLLs   += { "GRB" };
//==================================================
\endverbatim

<br>
At this point GRB is available from an external application (as Gleam). 
To have the item "GRB" in the sources menu of the GUI one can just edit 
the xml file containing the source definition adding the following lines:  
<br>
\verbatim 
<source name="GRB" >
<spectrum> <SpectrumClass name="GRBSpectrum"/> <use_spectrum/> </spectrum>
</source>
\endverbatim
<hr>


The /src/test/jobOptions.txt file holds information used for the
implementation of GRB algorithm. It doesn't contain any information 
regarding the physics of the GRB (that are all included in the GRBParam 
file), but it manages some options available for the GRB Algorithm.

\section jobOptions jobOptions files
There are two different jobOptions.txt file:

- <a href="../../src/test/test_jobOptions.txt>test_jobOptions.txt</a> is meant to be used in nightly 
builds together with test_GRB.exe

- <a href="../../src/test/jobOptions.txt>jobOptions.txt</a>:
This file is used to choose (by picking up/commenting out '\c #include' statements) between the 3 files: 
 GRBtestAlgOptions.txt, TDSreadFluxOptions.txt, and LatGRBAlertOptions.txt

- <a href="../../src/test/GRBtestAlgOptions.txt>GRBtestAlgOptions.txt</a>:
\param GRBTestAlg.source_name
passes the name of the GRB source, to be chosen among the ones defined in 
<a href="../../xml/GRB_user_library.xml>GRB_user_library.xml</a>
\param GRBTestAlg.background_name
passes the name of the background source, to be added on top of the GRB signal. 
It can be any spectrum defined in FluxSvc xml files.
\param GRBTestAlg.EvtMax
Maximum number of photon generated. Default is 100000
\param GRBTestAlg.savefile
"root" saves data in ROOT format, "ascii" in ASCII text file. 
Saving in both is also possible

- <a href="../../src/test/TDSreadFluxOptions.txt>TDSreadFluxOptions.txt</a>:
\param TDSreadFluxAlg.savefile
"root" saves data in ROOT format, "ascii" in ASCII text file. 
Saving in both is also possible

- <a href="../../src/test/LatGRBAlertOptions.txt>LatGRBAlertOptions.txt</a>:
\param LatGRBAlertAlg.nbckoff
Region threshold; determines when to start testing for false triggers. Default is 5.
\param LatGRBAlertAlg.mix
A value of 0 indicates that background mix has already been generated in file named by mixedFile field. Default is 0.
\param LatGRBAlertAlg.grbFile
Name of file listing events data.
\param LatGRBAlertAlg.backgroundFile
Name of file containing background data.
\param LatGRBAlertAlg.grbOffsetTime
Value to be used to offset events times. Default is 0.
\param LatGRBAlertAlg.mixedFile
Name of file containing background mixed data.

<hr>
\section requirements CMT requirements
\include cmt/requirements
<hr>
@section notes Release Notes
@include release.notes

<hr>
 */


