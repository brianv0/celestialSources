#include "DarkMatter/Profile.h"
#include "globals.h"


// compute densitiy (in GeV/cm^3) at 'rad' kpc  from the halo centre :
float Profile::evalDensity(float rad)
{
  float rc_cut = 1.e-5; //in kpc
  float r=(rad<rc_cut)?rc_cut:rad;
//   if (rad<rc_cut) 
//     r=rc_cut;
//   else 
//     r = rad;
  float esp = (m_beta-m_rgamma)/m_alpha;
  float rho = m_rho0 * pow(r/m_r0,-m_rgamma) * pow(1.+pow(m_r0/m_a,m_alpha),esp) / pow(1.+pow(r/m_a,m_alpha),esp);

  return rho;
}

float Profile::getRadius(float l,float l_int,float b_int)
{
  double l_halo = m_center.l();
  double b_halo = m_center.b();

  float R_halo=m_r0;
  float theta_int  = (90.-b_int)/180.*pi;
  float   phi_int  = (l_int-180.)/180.*pi;
  float sine_theta = sin(theta_int);
  float cosi_theta = cos(theta_int);
  float sine_phi   = sin(phi_int);
  float cosi_phi   = cos(phi_int);

  float  theta_halo = (90.-b_halo)/180.*pi;
  float    phi_halo = (l_halo-180.)/180.*pi;
  // Given position of the Halo (R_halo,theta_halo,phi_halo)
  float x_halo = R_halo*sin(theta_halo)*cos(phi_halo);
  float y_halo = R_halo*sin(theta_halo)*sin(phi_halo);
  float z_halo = R_halo*cos(theta_halo);

  float X   = x_halo-l*sine_theta*cosi_phi;
  float Y   = y_halo-l*sine_theta*sine_phi;
  float Z   = z_halo-l*cosi_theta;
  float rad = pow(X*X+Y*Y+Z*Z,0.5);
  return rad;
}
/////////////////////////////////////////////////////////
