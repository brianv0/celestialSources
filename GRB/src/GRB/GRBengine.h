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
  double getDurationFromBATSE(char* burst_type="Both");
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
