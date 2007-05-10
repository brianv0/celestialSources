/*! \mainpage microQuasar Package: simulating x-ray binary/microquasar source for GLAST
 * 
 *
 * \author Richard Dubois  richard@slac.stanford.edu 
 *
* @section introduction Introduction
*
* microQuasar is a package that simulates gamma-ray emission from x-ray binaries/microquasars.
* The model simulates three aspects:
* - sinusoidal modulation of the flux across the binary orbit period
* - two regions of different spectral index across the orbit
* - a cyclic "on" period for the high energy jet during the disk cycle
*

* A description for how the sinusoidal modulation works can be found at:
* <br>
* <a href="http://d0.phys.washington.edu/~burnett/glast/generate_periodic/oscilations.htm">http://d0.phys.washington.edu/~burnett/glast/generate_periodic/oscilations.htm </a>
*
* There are 16 parameters (case insensitive) to drive this model, supplied in the source xml file:
*
* flux - flux (p/s/cm^2)
* EMin - E(min) (Mev)
* EMax - E(max) (MeV)
* OrbitalPeriod - Orbital period (days)
* OrbitalModulation - orbital modulation (0 to 1)
* OrbitalPhase - orbital phase (0 to 1)
* SpectralOrbitalRegion1 - spectral index in orbital region 1 (positive)
* SpectralOrbitalRegion2(*) - spectral index in orbital region 2 (positive)
* OrbitalPhaseRegion1 - min orbital phase region 1 (0 to 1)
* OrbitalPhaseRegion2(*) - max orbital phase region 2 (0 to 1)
* DiskCycleDuration - disk cycle duration (days)
* DiskCycleFluctuation - disk cycle duration fluctuation (fraction: 0-1) - not implemented yet
* JetOnCycle(**) - start jet-on period in disk cycle (fraction: 0-1)
* JetOnCycleFluctuation(**) - fluctuation on jet start (fraction of jet-on start time: 0-1) - not implemented yet
* JetOnDuration(**) - duration of jet-on period (fraction disk cycle: 0-1)
* JetOnDurationFluctuation(**) - fluctuatian on jet-on duration (fraction of jet-on time: 0-1) - not implemented yet
*
*  (*) - optional if SpectralOrbitalRegion = -1; then only one region used
*  (**) - optional if JetOnCycle = 1; then no outbursts, just steady emission
*/

