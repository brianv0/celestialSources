 /*! \mainpage GRB Package: simulating a transient source
   
 This package contains three programs related to Gamma-Ray Bursts. The full physical simulator of GRB, based on the fireball model in the internal shocks configuration.

The phenomenological simulator of GRB signal. And the algorithm describing the LAT alert.

 GRB is now defined as a Gaudi algorithm, and it can be dynamically loaded within FluxSvc. GRBSpectrum is an example of rapid transient flux, and can be used for studyng the capability of Glast in the observation of rapid transient signals.
 
 \section physics Very short introduction to the physical model
<br>
The physical model we adopted is based on the \em fireball model of Gamma Ray Burst, for which a series of shells is injected in the circum burst medium with different Lorentz factor. When a faster shell reach a slower one a shock occours, and an accelerated electron distribution is obtained due to the shock acceleration process. Some of the energy dissipated during the shock is converted into a randomly oriented magnetic field. The electrons can loose their energy via synchrotron emission. The caratheristic synchrotron spectrum is boosted (and bemed) thancks to the Lorentz factor of the emittin matherial. The higher energy part of a GRB spectrum can be obtained keeping into account the possibility of Compton scattering of the synchrotron photons against the electron accellerated by the shock (Inverse Compton Scattering).
<br>
We adopted this model because it is able to well reproduce the fast time variability observed in the GRB signal. 
<br>
\section test How-to use the test programs

The GRB simulator can be use with several test program:

- test_GRB.exe
This executable tests the GRB algorithm. It initializes the GRB simulation,and it extracts photons according with the computed spectrum. To launch it type "test_GRB.exe". All the option are contained in the joboptions.txt file.

- test_GRBROOT.exe
This test program makes use of the \em ROOT graphical environment to display some plots regarding the simulated GRB. It shows the evolution of the flux with respect to the time, and it plots the integrated spectrum and the light curves.
to execute it type "test_GRBROOT.exe".
<br>

\section Howto How-to use the GRB spectrum in Gleam
<br>
GRB is an independent package, it is an external service and it has to be 
declared in the joboption file of FluxSvc. For example here there is the 
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
The file /src/test/GRBParam.txt holds information used for the description of the GRB. It can be changed without recompiling the entire package.
<br>
It MUST be in /src/test directory of the GRB package. 
<br>
The File GRBParam.txt contains information about:

-# Number of shells emitted by the engine (nshell).
-# The redshift of the source.
-# the total energy available in erg.
-# The initial separation from one shell to the subsequent.
-# The initial thickness of the shells.
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

\include jobOptions.txt

\section requirements requirements
\include requirements


 */


