 /*! \mainpage GRB Package: simulating a transient source
   
The GRB package has been developed to provide Gleam with a full fledged simulator of a physical transient source.
As such it has been designed to correctly interface FluxSvc, by providing spectrum classes derived from ISpectrum.
This package is also developed to study the physics of bursts and the capability of Glast in the observation of rapid 
transient signals. It is interfaced to ROOT for possible visualization.
 
 This package contains three parts related to Gamma-Ray Bursts:
 -# The full physical simulator of GRB, based on the fireball model in the internal shocks configuration.
 -# The phenomenological simulator of GRB signal. 
 -# The algorithm describing the LAT alert.

 
 \section physics Very short introduction to the physical model
<br>
We adopted a physical model based on the \em fireball model of Gamma Ray Burst, because it is able to well reproduce the fast time variability observed in the GRB signal: a series of shells is injected in the circum burst medium with different Lorentz factor. When a faster shell reach a slower one a shock occurs, and an accelerated electron distribution is obtained due to the shock acceleration process. Some of the energy dissipated during the shock is converted into a randomly oriented magnetic field. The electrons can loose their energy via synchrotron emission. The charateristic synchrotron spectrum is boosted (and beamed) thanks to the Lorentz factor of the emitting matherial. The higher energy part of a GRB spectrum can be obtained keeping into account the possibility of Compton scattering of the synchrotron photons against the electron accelerated by the shock (Inverse Compton Scattering).
<br>

\section phomenological Very short introduction to the phenomenological model

\section alert Very short introdction to the alert algorithm

<br>
\section test How-to use the test programs


The GRB simulator can be use with several test program:

- test_GRB.exe
This executable tests the GRB algorithm. It initializes the GRB simulation,and extracts photons according to the computed spectrum. To launch it type "test_GRB.exe". All the options are contained in the joboptions.txt file.

- test_GRBROOT.exe
This test program makes use of the \em ROOT graphical environment to display some plots regarding the simulated GRB. It shows the evolution of the flux with respect to time, and plots the integrated spectrum and the light curves.
to execute it type "test_GRBROOT.exe".
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


\section parameters GRB Parameters
<br>
The file /src/test/GRBParam.txt holds information used for the description of the GRB. 
It can be changed without recompiling the entire package.
<br>
It MUST be in /src/test directory of the GRB package. 
<br>
The File GRBParam.txt contains information about:

-# the number of shells emitted by the engine (nshell).
-# the redshift of the source.
-# the total energy available in erg.
-# the initial separation from one shell to the subsequent.
-# the initial thickness of the shells.
-# the minimum Lorentz factor of the shells
-# the maximum Lorentz factor of the shells

\verbatim
10		        // Nshell
1.0		        // Redshift
1.0e+52              // Etot
1.0e+8		// Radius of the shells [cm]
1.0e+8      	        // Thickness of the shells[cm]
100		        // Minumum Lorenz Factor
1000		        // Maximum Lorenz Factor
\endverbatim

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
\param same as above, replacing \c GRBTestAlg with \c FluxAlg.
- <a href="../../src/test/LatGRBAlertOptions.txt>LatGRBAlertOptions.txt</a>:
\param same as above, with \c LatGRBAlertAlg as the driving algorithm, plus
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
\section notes Release notes
doc/release.notes
<hr>


 */


