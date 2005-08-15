//////////////////////////////////////////////////
// File PulsarSim.h
// Header file for PulsarSim class
//////////////////////////////////////////////////

#ifndef PulsarSIM_H
#define PulsarSIM_H 

#include <vector>
#include <fstream>
#include <iomanip>
#include <ctime>
#include "TFile.h"
#include "TF1.h"
#include <TTree.h>
#include "PulsarConstants.h"
#include "SpectObj/SpectObj.h"
#include "facilities/Util.h"


/*! 
 * \class PulsarSim
 * \brief Class that contains the generation of the TH2D ROOT Histogram for Pulsar flux according to a selected model.
 *  
 * \author Nicola Omodei        nicola.omodei@pi.infn.it 
 * \author Massimiliano Razzano massimiliano.razzano@pi.infn.it
 *
 * This class creates the TH2D ROOT histogram that contains the differential photon flux (dN/dE/dt/dA) of the simulated
 * pulsar espressed in ph/keV/s/m2.
 * The user can specify the emission model. Now the default model is a phenomenological one based on observations of
 * known gamma-ray pulsars.
*/

class PulsarSim 
{
 public:
  
  //! Constructor of PulsarSim
  PulsarSim(std::string name, int seed, double flux, double enphmin, double enphmax, double period);

  //! Destructor of PulsarSim
  ~PulsarSim()
    {
      delete m_Nv;
    }
  
  //! Method that creates the TH2D histogram according to the phenomenological model.
  TH2D* PSRPhenom(double par0, double par1, double par2, double par3, double par4);

  // Returns a TH2D ROOT matrix that contains in every bin Nv*dE*dT*Aeff
  TH2D *PulsarSim::Nph(const TH2D *Nv);

  //! Returns the period of the pulsar 
  inline double Period(){return m_period;}

  //! Save a ROOT file with the TH2D ROOT histogram
  void SaveNv(TH2D *Nv);

  //! Save a TXT file with the Pulsar time profile
  void SaveTimeProfile(TH2D *Nv);
  
 private:
  
  //! Gathers all relevant constants for the simulation 
  double m_period;
  double m_flux;
  int m_Tbin;
  int m_numpeaks;
  double m_enphmin, m_enphmax;
  int m_seed;
  std::string m_name;
  TH2D *m_Nv;
};

#endif


