/*!
  \class GRBobsSim
  \brief Simulator engine of a GRB source.
 
  \author Nicola Omodei       nicola.omodei@pi.infn.it 
 
*/

#include "GRBobsConstants.h"
#include "GRBobsPulse.h"
#include "GRBobsengine.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#ifndef GRBobsSIM_H
#define GRBobsSIM_H 1

class GRBobsSim 
{
 public:

  GRBobsSim(GRBobsParameters *params);
  //! destructor
  ~GRBobsSim()
    {
      delete m_GRBengine;
      delete m_Nv;
    }
  
  void GetUniqueName(const void *ptr, std::string & name);
    
  /*!
   * \brief Starts the GRBobs simulation
   *
   * Initialize the simulation
   */
  TH2D* MakeGRB();

  /*! Compute the Flux, as a function of time. It returns a matrix.
   * \param time is the time in which the spectrum is calculated.
   */
  TH2D *GRBobsSim::Nph(const TH2D *Nv);
  inline std::pair<double,double> GRBdir(){return m_GRBengine->GetDirection();}
  inline double Tmax(){return m_tfinal;}
  void SaveNv();
  void GetGBMFlux();
 private:
  
  //! Gathers all relevant constants for the simulation 
  GRBobsParameters *m_params;
  GRBobsengine  *m_GRBengine;
  double m_tfinal;
  double m_fluence;
  int m_tbin;
  TH2D *m_Nv;
};

#endif


