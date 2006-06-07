#ifndef PROFILE_H
#define PROFILE_H

#include "astro/SkyDir.h"

class Profile 
{
public:
  //void Profile(const astro::SkyDir& bincenter, float distance, const char* type, const float rho0);
  Profile(const astro::SkyDir& center, float a, float alpha, float beta, float rgamma, float rho0,float r0)
  {
    m_center = center;
    m_a = a;
    m_alpha = alpha;
    m_beta =  beta;
    m_rgamma =  rgamma;
    m_rho0 =  rho0;
    m_r0 =  r0;
  }
  float evalDensity(float);
  float getRadius(float l,float l_int,float b_int);

  //private:
  astro::SkyDir m_center;
  float distance;
  float m_a;
  float m_alpha;
  float m_beta;
  float m_rgamma;
  float m_r0;
  float m_rho0;
};

#endif
