#ifndef GRBmanager_H
#define GRBmanager_H

#include <vector>
#include <string>
#include <map>
#include <cmath>
#include "FluxSvc/ISpectrum.h"
#include "facilities/Observer.h"
#include "src/GPS.h"
#include "CLHEP/Random/RandomEngine.h"
#include "GRBSpectrum.h"

class ISpectrum;

class GRBmanager : public ISpectrum
{

 public:
  
  //! Constructor: takes a file with some parameters as argument.
  /*! \param params file of parameters 
   */
  GRBmanager(const std::string& params);
  /*! Destructor
   */
  virtual  ~GRBmanager(); 
  
  double flux(double time)const; //{return m_GRB->flux(time);}
  double interval(double time);//{return m_GRB->interval(time);}//const;
  double solidAngle() const;//{return m_GRB->solidAngle();}
  
  inline std::pair<double,double> 
    dir(double energy, HepRandomEngine* engine){return m_GRB->dir(energy, engine);} 


  float operator() (float u) const;
  double energySrc(HepRandomEngine*, double /*time*/ );
  //! inherited from Spectrum
  std::string title() const {return "GRBmanager";} 
  const char * particleName() const {return "gamma";}
  const char * nameOf() const {return "GRBmanager";}
  float parseParamList(std::string input, int index);  
 private:
  
  GRBSpectrum* m_GRB;
  const std::string& m_params;
  float m_FirstTime;
  float m_initialTime;
  float m_timeToWait;
  float m_endTime;
};
#endif
