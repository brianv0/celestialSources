/*! \mainpage PulsarSpectrum Package: simulating a Pulsar transient source for GLAST
 * 
 *
 * \author Nicola Omodei        nicola.omodei@pi.infn.it 
 * \author Massimiliano Razzano massimiliano.razzano@pi.infn.it
 *
 *PulsarSpectrum is a package that simulates gamma-ray emission from pulsars within different models and scenarios. It is based on the same structure of the GRB package contained in celestialSources. It is fully compatible with the GLAST LAT software, and can be easily interfaced with Gleam or ObservationSim in order to produce photons from gamma ray pulsars for studies on science performances and data analysis.
 * The Simulator is also able to compute the timing decorrections due to period change and barycentri decorections.
 *As the GRB simulator, PulsarSpectrum creates a TH2D ROOT histogram that contains the flux of the simulated sources (ph/m2/kev/s) vs. energy (keV) and time (s.). The time profile of the source is a 1 (or 2) peaked curve (random generated), and the spectral profile is created according to the emission model choosen by the user. This TH2D histogram is then available to Gleam and ObservationSim simply through the flux package.
The source code for PulsarSpectrum is included in the celestialSources package under the directory named Pulsar. The test program is named test_PulsarROOT.exe and produce some plots of an example of pulsar whose parameters are set in the code of the executable.

* At the moment there's included one default emission model, based on a phenomenological analytical form for the spectrum, based on Nel and DeJager (1995, see Ref.) and of DeJager (2002, see Ref.) that takes into account the spectral high energy cutoff. Using this model, in this tutorial is shown an example that describes how to create pulsar sources with a Polar Cap or Outer Gap scenario.

* A brief tutorial on the use of PulsarSpectrum can be found at:
* <br>
* <a href="#dejager02">http://www.pi.infn.it/~razzano/Pulsar/PulsarSpTutor/PulsarSpTutor.htm </a>
*/
