#include "globals.h"
#include "DarkMatter/Profile.h"
#include "astro/SkyFunction.h"

class LineOfSight : virtual public astro::SkyFunction
{
public:
  LineOfSight(Profile& prof);
  ~LineOfSight() {delete m_profile;}
  double operator()(const astro::SkyDir& bincenter) const;
  //void setIntegrationSteps();
  float integral1(float l_obs,float b_obs,float time);
  void setBinSize(const float dl,const float db)
    {
      m_dl=dl;
      m_db=db;
      float solangle=2*pi*(1-cos(db*pi/180.));
      r_inside=pow(solangle/pi,0.5)*m_profile->m_r0;
      std::cout<<"r_inside="<<r_inside<<" kpc"<<std::endl;
      step1   = r_inside/100.;
    }

private:
  Profile* m_profile;
  float step0;
  float step1;
  float step2;
  float step3;
  float r_step2;
  float r_step3;
  float r_inside;
  float j_inside;

  float m_dl;
  float m_db;
};

