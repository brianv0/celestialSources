#include "LineOfSight.h"

LineOfSight::LineOfSight(Profile& prof)
{
  m_profile = &prof;
  float R_halo = m_profile->m_r0;

  r_inside=0.1; //this value is reevaluated by setBinSize()

  // def of steps for in tegration :
  step1   = r_inside/100.;
  step2   = 0.01;
  r_step2 = 0.1;
  step3   = 0.1;
  r_step3 = 1.;

  if (step1>step2)
    std::cout<< "Problem : step1 > step2"<<std::endl;

  if (R_halo>r_step2)
    {
      if (R_halo>r_step3)
		 step0=step3;
      else
	step0=step2;
    }	
  else
    step0=step1;

  float j_inside=0; //looks useless
  
}

float LineOfSight::integral1(float l_obs,float b_obs,float time)
{

  float delta_l=m_dl;
  float delta_b=m_db;

  float R_halo=m_profile->m_r0;
    
  float numX1=time+1;
  float numY1=time+1;
  //jmap1=num.zeros((numX1,numY1),type=num.Float64);
  float b_int=b_obs-delta_b/2.;
  int ix1=0;

  float integral=0.;
  while (b_int<=b_obs+delta_b/2.)
    {
      float l_int=l_obs-delta_l/2.;
      int iy1=0;
      while (l_int<=l_obs+delta_l/2.)
	{
	  float step=step0;
	  float l_coor=0.;
	  float j=0.;
	  int flag=0; //Show if the line of sight go through the inside region
	  while (l_coor<R_halo+10)
	    {
	      float l=l_coor+step*0.5;
	      float rad=m_profile->getRadius(l,l_int,b_int);
	      float rho=m_profile->evalDensity(rad);

	      if (rad>r_step3)   // r_step3 = 1    kpc
		step=step3;    //   step  = 0.1  kpc
	      else
		{
		  if (rad>r_step2) // r_step2 = 0.1  kpc
		    step=step2;    //   step  = 0.01 kpc
		  else
		    step=step1;    //   step  = r_inside/100.
		}	      
	      j+=step*pow(rho,2);
	      l_coor+=step;
	    }
	  integral+=j; 
	  //jmap1[ix1,iy1]=j;
	  l_int+=delta_l/time;
	  iy1+=1;
	}
      b_int+=delta_b/time;
      ix1+=1;
    }
  return integral;
  //return integral*delta_l/time*delta_b/time*pow(deg2rad,2)*kpc2cm;
}

double LineOfSight::operator()(const astro::SkyDir& bincenter)const
{

   float l_obs=bincenter.l();
   float b_obs=bincenter.b();
   float delta_l=m_dl;
   float delta_b=m_db;

   float R_halo=m_profile->m_r0;
   float l_halo=m_profile->m_center.l();
   float b_halo=m_profile->m_center.b();

   float integral=0.;
   float sig19=0.;

  if ((fabs(l_obs-l_halo)<delta_l)&&(fabs(b_obs-b_halo)<delta_b))
    {
      std::cout<<"computing center in better grid"<<std::endl;
      int time=20;
      float rint=integral1(l_obs,b_obs,time);
      float rDeltaOmega=delta_l/time*delta_b/time*pow(deg2rad,2)*kpc2cm;
      integral+=rint;
      sig19   +=rint*rDeltaOmega;
      //jmap[ix][iy]=rint;
      std::cout<< "central value of Sigma 19= "<<b_obs<<" "<<l_obs<<" "<<sig19/1.e+19<<std::endl;
    }
  else
    {
      float step=step0;
      float l_coor=0.;
      float j=0.;
      int flag=0; //Show if the line of sight go through the inside region
      while (l_coor<R_halo+10.) 
	{
	  float l=l_coor+step*0.5;
	  float rad=m_profile->getRadius(l,l_obs,b_obs);
	  float rho=m_profile->evalDensity(rad);
	  if (rad>r_step3)     // r_step3 = 1    kpc
	    step=step3;        //   step  = 0.1  kpc
	  else
	    {
	      if (rad>r_step2) // r_step2 = 0.1  kpc
		step=step2;    //   step  = 0.01 kpc
	      else
		step=step1;    //   step  = r_inside/100.
	    }
	  
	  if (rad<r_inside)
	    {
	      flag+=1;
	      l_coor+=step*10.;
	    }
	  else
	    {
	      j+=step*pow(rho,2);
	      l_coor+=step;
	    }
	}
      if (flag>0)
	{
	  j+=j_inside;
	  std::cout<< "wrong"<<std::endl;
	}
      // in fact there is a factor pow(R_halo,2)/R_halo/R_halo :
      float rDeltaOmega=delta_l*delta_b*pow(deg2rad,2)*kpc2cm;
      integral   +=j;
      sig19      +=j*rDeltaOmega;
    } 
  return integral;
}
