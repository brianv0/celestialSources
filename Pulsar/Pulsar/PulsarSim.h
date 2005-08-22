#include "PulsarConstants.h"
#include "SpectObj/SpectObj.h"
#include <vector>

#ifndef PulsarSIM_H
#define PulsarSIM_H 

class PulsarSim 
{
 public:
  
  //! The simulation can be initialized by setting the seed of the Random engine generator. If seed is 0, a Random Pulsar is computed.
  PulsarSim(double fluence,double period, int numpeaks);
  //! destructor
  ~PulsarSim()
    {
      delete m_Nv;
    }
  
  /*!
   * \brief Starts the Pulsar simulation
   *
   * Initialize the simulation
   */
  TH2D* PSRPhenom(double par1, double par2, double par3, double par4);

  /*! Compute the Flux, as a function of time. It returns a matrix.
   * \param time is the time in which the spectrum is calculated.
   */
  TH2D *PulsarSim::Nph(const TH2D *Nv);
  inline double Period(){return m_period;}
  void SaveNv(TH2D *Nv);
 private:
  
  //! Gathers all relevant constants for the simulation 
  double m_period;
  double m_fluence;
  int m_numpeaks;
  TH2D *m_Nv;
};

#endif


