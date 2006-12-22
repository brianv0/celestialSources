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
* There are 16 parameters to drive this model, supplied in the source xml file:
*
* - 0 - flux (p/s/cm^2)
* - 1 - E(min) (Mev)
* - 2 - E(max) (MeV)
* - 3 - Orbital period (days)
* - 4 - orbital modulation (0 to 1)
* - 5 - orbital phase (0 to 1)
* - 6 - spectral index in orbital region 1 (positive)
* - 7 - spectral index in orbital region 2 (positive)
* - 8 - min orbital phase region 1 (0 to 1)
* - 9 - max orbital phase region 1 (0 to 1)
* - 10 - disk cycle duration (days)
* - 11 - disk cycle duration fluctuation (fraction: 0-1) - not implemented yet
* - 12 - start jet-on period in disk cycle (fraction: 0-1)
* - 13 - fluctuation on jet start (fraction of jet-on start time: 0-1) - not implemented yet
* - 14 - duration of jet-on period (fraction disk cycle: 0-1)
* - 15 - fluctuatian on jet-on duration (fraction of jet-on time: 0-1) - not implemented yet
*
*/

