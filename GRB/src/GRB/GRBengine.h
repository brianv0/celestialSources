/*!
 * \class GRBengine
 *
 * \brief This class permits to generate different kind of GRBs.
 * 
 * \author Nicola Omodei       nicola.omodei@pi.infn.it 
 * \author Johann Cohen-Tanugi johann.cohen@pi.infn.it
 *
 */

#ifndef GRBENGINE_H
#define GRBENGINE_H 1

#include "GRBConstants.h"
#include "GRBShock.h"

class GRBengine
{
  
 public:
  GRBengine(Parameters *params);
  
  //  void   EraseShocksVector();
  
  std::vector<GRBShock*> CreateShocksVector(const int N,
					    const double L,
					    const double l,
					    const double etot
					    );
  double GenerateGamma(double gammamax,double gammamin);
  //double getDistance();
  inline std::pair<double,double> GetDirection(){return m_dir;}
 private:
  std::pair<double,double> m_dir;
  Parameters *m_params;
};

#endif
