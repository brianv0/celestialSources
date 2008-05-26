/*! \mainpage PulsarSpectrum Package: simulating a Pulsar transient source for GLAST
 * 
 *
 * \author Nicola Omodei        nicola.omodei@pi.infn.it 
 * \author Massimiliano Razzano massimiliano.razzano@pi.infn.it
 *
 *PulsarSpectrum is a package that simulates gamma-ray emission from pulsars within different models and scenarios. It is based on the same structure of the GRB package contained in celestialSources. It is fully compatible with the GLAST LAT software, and can be easily interfaced with Gleam or ObservationSim in order to produce photons from gamma ray pulsars for studies on science performances and data analysis.
 * The Simulator is also able to compute the timing decorrections due to period change and barycentric decorections.
 *As the GRB simulator, PulsarSpectrum creates a TH2D ROOT histogram that contains the flux of the simulated sources (ph/m2/kev/s) vs. energy (keV) and time (s.). The time profile of the source is a 1 (or 2) peaked curve (random generated), and the spectral profile is created according to the emission model choosen by the user. This TH2D histogram is then available to Gleam and ObservationSim simply through the flux package.
The source code for PulsarSpectrum is included in the celestialSources package under the directory named Pulsar. The test program is named test_PulsarROOT.exe and produce some plots of an example of pulsar whose parameters are set in the code of the executable.
*
* 2 simulation models are included:
* - PSRPhenom, based on an analytical formula for the spectrum based on Nel and DeJager (1995 ;
* - PSRShape, that allows the user to simulate arbitrary phase-energy models;
*
* Environment variables;
* In order to customize PulsarSpectrum it is possible to set some environment variables, here described:
*
* PULSAR_OUTPUT_LEVEL:
* Sets the level of output log
*   0 - No log
*   1 - Only pulsar model information
*   2 - Pulsar model information and warnings
*
* PULSAR_OUT_NV
* If not defined, no Nv ROOT histogram is saved, otherwise if is set (to whatever value) the Nv ROOT file is written 
*   
* PULSAR_OUT_BARY
* If not defined, no file with barycentric corrections is saved, otherwise if is set (to whatever value) a txt file called PulsarNameBaryCorr.txt is created for each pulsar called PulsarName;
*
* PULSAR_OUT_BIN
* If not defined, no file with binary demodulation is saved, otherwise if is set (to whatever value) a txt file called PulsarNameDemodBin.txt is created for each pulsar called PulsarName;
*
* PULSAR_OUT_TNOISE
* If not defined, no file with timing noise residual is saved, otherwise if is set (to whatever value) a txt file called PulsarNameTNoiseLog.txt is created for each pulsar called PulsarName;
*
* PULSAR_NO_DB
* If not defined the txt files containing the database are set, otherwise if this env variable is set (to whatever value) no output .txt database file will be written 
*
* PULSAR_EPH
* Ephemerides label in the output D4 of the simulations. If not defined, the defaults ephemerides used is DE405
*
*
* A brief tutorial on the use of PulsarSpectrum can be found at:
* <br>
* <a href="#dejager02">http://www.pi.infn.it/~razzano/Pulsar/PulsarSpTutor/PulsarSpTutor.htm </a>
*/
