/*!
 * \class GRBengine
 *
 * \brief This class permits to generate different kind of GRBs.
 * 
 * There are different type of engine:
 *
 * - Type 0: GRBengine uses some relations between the physical
 *    parameters, such as the dimension of the shocked material,
 *    and some quantities, such as the duration of a spike,
 *    and position of the break energy in the flux.
 *    In particular, if \f$t_r\f$ and \f$t_d\f$ are the rise and 
 *    the decay time of a single spike, \f$\delta_r\f$ is the thickness 
 *    of a shell, \f$E_i\f$ is the internal energy available in a shock,
 *    then the estimated relations are: 
 *
 *    \f$t_r\propto\displaystyle{\delta_r\over \Gamma c}\f$.
 *
 *    \f$t_d\propto\displaystyle{jet_r^2 \delta_r \over \Gamma^2 E_i  }\f$.
 *
 *    \f$E_b\propto\displaystyle{\Gamma^3 \sqrt{E_i/\delta_r}\over jet_r  }\f$.
 *
 *    From this equations it is easy to estimate three parameters, 
 *    living free the fourth (\f$E_i\f$).
 *
 * - Type 1: GRBengine will create a Shocks sequence, 
 *    according with the physical parameter chosen in the configuration 
 *    file The file  <a href="../../src/test/GRBParam.txt"> GRBParam.txt </a>,
 *    and stored in GRBConstants.
 *    GRBengine will read the number of shocks,
 *    the duration of the burst,
 *    the energy, the Lorentz factor, and the thickness of the shell.
 *    Depending on the Shell_type, GRBengine will read the radios of the shell,
 *    or the jet angle and the jet radius.

 * - Type 2: GRBengine will first create two Shells for each shock,
 *    retrieving the characteristic parameters of GRBShells from GRBConstants,
 *    and then will compute their collisions giving up to the GRBShock 
 *    sequence.
 
 * \e GRBengine \e does \e the \e following:
 * - Reads the parameter file.
 * - Depending on the particular type of engine, fills a shocks vector
 * (that is a vector of GRBShock)
 * - Sets the distance and the position in the sky of the GRB.
 * 
 * \author Nicola Omodei       nicola.omodei@pi.infn.it 
 * \author Johann Cohen-Tanugi johann.cohen@pi.infn.it
 *
 */


#include "GRBShock.h"
#include "GRBConstants.h"
#include <vector>

#ifndef GRBENGINE_H
#define GRBENGINE_H 1

class GRBengine
{
 public:

  GRBengine(GRBConstants *myParam);
  //GRBengine(){;}
  ~GRBengine(){;}
  /*!The observed duration from CGRO/BATSE has been fitted with a double power law.
   * This method returns a number chosen from this distribution.
   * \param burst_type is the type of GRB: "Short" only short bursts 
   * are considered, if it is "Long", only long bursts. 
   * "Both" means that the two distribution
   * are considered, selecting 30% short bursts and 70% long bursts.
   */
  double getDurationFromBATSE(char* burst_type="Both");
  double GetTime(double duration,int nshok,double rnd);
  inline double getDuration(){return m_duration;}
  inline double getDistance(){return m_distance;}
  inline std::pair<double,double> getDirection(){return m_direction;}
  std::vector<GRBShock> getShocksVector(){return theShocks;}
 private:
  double m_duration;
  double m_distance;
  std::pair<double,double> m_direction;
  std::vector<GRBShock> theShocks;
  
};
#endif
