#include <iostream>
#include <math.h>
#include <sstream>
#include <string>
#include <stdexcept>

#include "tip/IFileSvc.h"
#include "tip/Image.h"

#include "fitsio.h"

#define kpc2cm  3.080e+21
#define pi      3.1415926
#define parsec  3.085677E16  /* m */
#define deg2rad pi/180.     

static std::string s_fitsRoutine("");

float r0  ; //in kpc
float rho0; //in GeV/cm^3
float a   ; //in kpc
        
float alpha;
float beta ;
float rgamma;

float rc_cut ; //in kpc

float dl;
float db;

float l_halo;
float b_halo;

float step0;
float step1;
float step2;
float step3;
float r_step2;
float r_step3;
float r_inside;
float j_inside;

float jmap[1000][1000];

float evalDensity(float);
float getRadius(float, float, float);
float integral1(float, float, float);
float LOSIntegral(float, float ,float ,float );



void fitsReportError(int status, std::string routine="") {
   if (status == 0) {
      return;
   }
   if (routine == "") {
      routine = "haloProfile::" + s_fitsRoutine;
   }
   fits_report_error(stderr, status);
   std::ostringstream message;
   message << routine << ": CFITSIO error " << status;
   throw std::runtime_error(message.str());
}

void create2DMap(std::string & filename)
{
  //using tip to create the fits file
  tip::IFileSvc::instance().createFile("new_image.fits");
  // Open new image for writing.
  int ndimX= int(360/dl);
  int ndimY=int(180/db);
  tip::ImageBase::PixelCoordinate img_dims;
  img_dims.push_back(ndimX);
  img_dims.push_back(ndimY);
  tip::IFileSvc::instance().appendImage("new_image.fits", "Primary", img_dims);
  
  tip::Image* tipImage = tip::IFileSvc::instance().editImage("new_image.fits", "Primary");
  tip::ImageBase::PixelCoordinate pixel;
  pixel.push_back(0);
  pixel.push_back(0);      
  tipImage->setPixel(pixel,0);

  for(int ix=0;ix<ndimX;ix++){
    for(int iy=0;iy<ndimY;iy++){
      //      tip::ImageBase::PixelCoordinate pixel;
      pixel.push_back(ix);
      pixel.push_back(iy);      
      tipImage->set(pixel,jmap[ix][iy]);
    }
  }



  //brutal cfitsio  call
  s_fitsRoutine="create2DMap";
  fitsfile * fptr = 0;
  char *file = const_cast<char *>(filename.c_str());
  int status = 0;
  
  fits_create_file(&fptr, file, &status);
  fitsReportError(status);
  
  long dims[2] = { 40,40 };
  fits_create_img(fptr, FLOAT_IMG, 2, dims, &status);
  fitsReportError(status);

  fits_close_file(fptr, &status);
  fitsReportError(status);

}

//----------------------------------------------------------------------
// define the profile
void Profile()
{
  r0   =  8.;  //in kpc
  rho0 =  0.3; //in GeV/cm^3
  a    = 20.;  //in kpc
        
  alpha  = 1.;
  beta   = 3.;
  rgamma = 1.;

  rc_cut = 1.e-5; //in kpc

  //Width of the bins Map:
  //dl = 0.1;
  //db = 0.1;
  dl = 0.01;
  db = 0.01;

  std::cout<<"Width of the bins Map: : "<<std::endl;
  std::cout<< "db = "<<db<<" dl = "<<dl<<std::endl;
 
  //halo position in Galactic Coordinate System
  l_halo = 0.;
  b_halo = 0.;

  float solangle=2*pi*(1-cos(db*pi/180.));

  float r_inside=pow(solangle/pi,0.5)*r0;
  std::cout<<"r_inside="<<r_inside<<" kpc"<<std::endl;

  // def of steps for in tegration :
  step1   = r_inside/100.;
  step2   = 0.01;
  r_step2 = 0.1;
  step3   = 0.1;
  r_step3 = 1.;

  if (step1>step2)
    std::cout<< "Problem : step1 > step2"<<std::endl;

  if (r0>r_step2)
    {
      if (r0>r_step3)
		 step0=step3;
      else
	step0=step2;
    }	
  else
    step0=step1;

  float j_inside=0; //looks useless
}
//----------------------------------------------------------------------
// compute densitiy (in GeV/cm^3) at 'rad' kpc  from the halo centre :
float evalDensity(float rad)
{
  float r;
  if (rad<rc_cut) 
    r=rc_cut;
  else 
    r = rad;
  float esp = (beta-rgamma)/alpha;
  float rho = rho0 * pow(r/r0,-rgamma) * pow(1.+pow(r0/a,alpha),esp) / pow(1.+pow(r/a,alpha),esp);

  return rho;
}
//----------------------------------------------------------------------
float getRadius(float l,float l_int,float b_int)
{
  float theta_int  = (90.-b_int)/180.*pi;
  float   phi_int  = (l_int-180.)/180.*pi;
  float sine_theta = sin(theta_int);
  float cosi_theta = cos(theta_int);
  float sine_phi   = sin(phi_int);
  float cosi_phi   = cos(phi_int);
  //----------------
  float  theta_halo = (90.-b_halo)/180.*pi;
  float    phi_halo = (l_halo-180.)/180.*pi;
  // Given position of the Halo (r0,theta_halo,phi_halo)
  float x_halo = r0*sin(theta_halo)*cos(phi_halo);
  float y_halo = r0*sin(theta_halo)*sin(phi_halo);
  float z_halo = r0*cos(theta_halo);
  //----------------
  float X   = x_halo-l*sine_theta*cosi_phi;
  float Y   = y_halo-l*sine_theta*sine_phi;
  float Z   = z_halo-l*cosi_theta;
  float rad = pow(X*X+Y*Y+Z*Z,0.5);
  return rad;
 }
//----------------------------------------------------------------------
float integral1(float l_obs,float b_obs,float time)
{
  float delta_l=dl;
  float delta_b=db;
    
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
	  while (l_coor<r0+10)
	    {
	      float l=l_coor+step*0.5;
	      float rad=getRadius(l,l_int,b_int);
	      float rho=evalDensity(rad);

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
//----------------------------------------------------------------------
float LOSIntegral(float lmin,float lmax,float bmin,float bmax)
{
  //  std::cout<< "Contour de la map : "<<std::endl;
  //  std::cout<< "lmin ="<<lmin<<" lmax="<<lmax<<std::endl;
  //  std::cout<< "bmin ="<<bmin<<" bmax="<<bmax<<std::endl;
  
  float delta_l=dl;
  float delta_b=db;
  
  int numX=int((bmax-bmin)/delta_b)+1;
  int numY=int((lmax-lmin)/delta_l)+1;

  for(int i=0;i<1000;i++)
    {
      for(int j=0;j<1000;j++)
	jmap[i][j]=0;
    }
  
  float b_obs=bmin+db/2.;
  int ix=0;
  float integral=0.;
  float sig19=0.;
  while (b_obs<bmax)
    {
      float l_obs=lmin+dl/2.;
      int iy=0;
      while (l_obs<lmax)
	{ 
	  if ((fabs(l_obs-l_halo)<delta_l)&&(fabs(b_obs-b_halo)<delta_b))
	    {
	      std::cout<<"computing center in better grid"<<std::endl;
              int time=20;
              float rint=integral1(l_obs,b_obs,time);
              float rDeltaOmega=delta_l/time*delta_b/time*pow(deg2rad,2)*kpc2cm;
	      integral+=rint;
	      sig19   +=rint*rDeltaOmega;
	      jmap[ix][iy]=rint;
	      //	      std::cout<< "central value of Sigma 19= "<<b_obs<<" "<<l_obs<<" "<<sig19/1.e+19<<std::endl;

	      l_obs+=delta_l;
	      iy+=1;
	      std::cout<< "central value of Sigma 19= "<<b_obs<<" "<<l_obs<<" "<<sig19/1.e+19<<std::endl;
	    }
	  else
	    {
	      float step=step0;
	      float l_coor=0.;
	      float j=0.;
	      int flag=0; //Show if the line of sight go through the inside region
	      while (l_coor<r0+10.) 
		{
		  float l=l_coor+step*0.5;
		  float rad=getRadius(l,l_obs,b_obs);
		  float rho=evalDensity(rad);

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

	      // in fact there is a factor pow(r0,2)/r0/r0 :
              float rDeltaOmega=delta_l*delta_b*pow(deg2rad,2)*kpc2cm;
	      integral   +=j;
	      sig19      +=j*rDeltaOmega;
	      jmap[ix][iy]=j;
	      l_obs+=delta_l;
	      
	      iy+=1;
	    }
	}
      b_obs+=delta_b;
      ix+=1;
    }
  
  //  std::cout<<std::endl<<"integral   = "<<integral/1.e+19<<std::endl<<std::endl;
  //  std::cout<<std::endl<<"Sigma19    = "<<sig19   /1.e+19<<std::endl<<std::endl;

  float rtest=0.;
  for(int i=0;i<1000;i++)
    {
      for(int j=0;j<1000;j++)
	rtest+=jmap[i][j];
    }

  //  std::cout<<"rtest      = "<<rtest/1.E19<<std::endl;

  return integral;
}

//----------------------------------------------------------------------
int main(void)
{
  Profile(); //profile definition

  float step=dl;
  int nbins=1000;
  float jmap=0.;
  for(int r=0;r<nbins;r++){
    float radius=r*step;
    jmap = LOSIntegral(radius-step/2.,radius+step/2.,radius-step/2.,radius+step/2.);
    std::cout<<radius<<"\t"<<jmap<<std::endl;    
  }
  //float jmap = LOSIntegral(-1.,1.,-1.,1.);
  //float jmap = LOSIntegral(-1.41,1.41,-1.41,1.41);
  //float jmap = LOSIntegral(-3.,3.,-3.,3.);
  //float jmap = LOSIntegral(-15.,15.,-15.,15.);
  //float jmap = LOSIntegral(-30.,30.,-30.,30.);

//   std::string filename("test.fits");
//   create2DMap(filename);

  std::cout<< "done"<<std::endl;
  return 0;
}
//----------------------------------------------------------------------
