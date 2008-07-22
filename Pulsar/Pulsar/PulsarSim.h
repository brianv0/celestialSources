/**
 * @file PulsarSim.h
 * @brief Class of the core of Pulsar Simulator, PulsarSim.cxx
 * @ author Massimiliano Razzano (massimiliano.razzano@pi.infn.it
 * @ author Nicola Omodei (nicola.omodei@pi.infn.it
 *
 * $Header$
 */
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
 * This class creates the TH2D ROOT histogram that contains the differential photon flux (dN/dE/dt/dA)
 * of the simulate pulsar espressed in ph/keV/s/m2.
 *
 * The user can specify the emission model. The two classes of model that can be used are:
 *
 * - PSRPhenom - A phenomenological model with an analitical spectrum (Nel & De Jager 1995);
 * - PSRShape  - A model that allow the user to use an arbitrary 2D ROOT spectrum model;
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

  //! Method that creates the TH2D histogram according to an external pulsar 2-d Shape
  TH2D* PSRShape(std::string ModelShapeName="PulsarShape", int NormalizeFlux=0);

  // Returns a TH2D ROOT matrix that contains in every bin Nv*dE*dT*Aeff
  TH2D *Nph(const TH2D *Nv);

  //! Returns the period of the pulsar 
  inline double Period(){return m_period;}

  //! Save a ROOT file with the TH2D ROOT histogram
  void SaveNv(TH2D *Nv);

  //! Save a TXT file with the Pulsar time profile
  void SaveTimeProfile(TH2D *Nv);

  //! Write to an output log file
  void WriteToLog(std::string Line);
  
  //!Flag for Output level
  int m_OutputLevel;

 private:

  //! Pulsar name
  std::string m_name;

  //! Pulsar Period
  double m_period;

  //! Pulsar flux
  double m_flux;

  //! Number of bin in time
  int m_Tbin;

  //! Number of peaks
  int m_numpeaks;

  //! Minimum and maximum energy of the extracted photon
  double m_enphmin, m_enphmax;

  //!Random seed
  int m_seed;

  //! PULSARDATA directory
  std::string m_pulsardata_dir;

  //!Outpur 2D ROOT histogram
  TH2D *m_Nv;

  //! output log filename
  std::string m_LogFileName;

};

#endif


