#include <iostream>
#include <math.h>
#include <sstream>
#include <string>
#include <stdexcept>

#include "map_tools/SkyImage.h"
#include "astro/SkyDir.h"

#include "DarkMatter/Profile.h"
#include "LineOfSight.h"

#include "fitsio.h"

static std::string s_fitsRoutine("");

void fitsReportError(int status, std::string routine="") {
   if (status == 0) {
      return;
   }
   if (routine == "") {
      routine = "makeMap::" + s_fitsRoutine;
   }
   fits_report_error(stderr, status);
   std::ostringstream message;
   message << routine << ": CFITSIO error " << status;
   throw std::runtime_error(message.str());
}

void testRadialDist(float step, int nbins, const LineOfSight los)
{
  for(int r=0;r<nbins;r++)
    {
      float radius = r*step;
      astro::SkyDir pixel_center=astro::SkyDir(radius,radius,astro::SkyDir::GALACTIC);
      std::cout<<radius<<"\t"<<los(pixel_center)<<std::endl;
    }
}

//----------------------------------------------------------------------
int main(void)
{
  ////////////////////////////////////////
  //SkyDir object defining the halo center
  float l_halo=0.;
  float b_halo=0.;
  astro::SkyDir halo_center = astro::SkyDir(l_halo,b_halo,astro::SkyDir::GALACTIC);

  ////////////////////////////////////////
  //Define a DM profile
  float r0     = 8.;  //in kpc
  float rho0   = 0.3; //in GeV/cm^3
  float a      = 20.;  //in kpc
  float alpha  = 1.;
  float beta   = 3.;
  float rgamma = 1.;
  Profile prof(halo_center,a,alpha,beta,rgamma,rho0,r0);

  ////////////////////////////////////////
  //Define the LineOfSight functor, derived from SkyFunction
  LineOfSight los(prof);
  //SkyImage only supports square pixels
  float pixel_size = 0.5;
  los.setBinSize(pixel_size, pixel_size);

  double fov=360.;
  std::string filename("output3.fits");

  // pixel_size [0.5] degree size of indivitual pixel
  // fov [20] (degrees) size of field of view, square if <90, full sky if>90
  // 1 = layers [1] number of layers to allocate (for energy scale)
  map_tools::SkyImage skyImage(halo_center,filename,pixel_size,fov,1,"CAR",true);
  skyImage.fill(los,0);

  std::cout<<"Sigma19: "<<los.getSigma19()<<std::endl;
  std::cout<<"Integral (without solid angle): "<<los.getIntegral()<<std::endl;

  s_fitsRoutine="create2DMap";
  fitsfile * fptr = 0;
  char *file = const_cast<char *>(filename.c_str());
  int status = 0;
  fits_open_file(&fptr, file, READWRITE,&status);
  fitsReportError(status);
  //fits_write_key(fptr, TFLOAT, "RZERO", &r0,"Distance to the sun in kpc", &status);
  fitsReportError(status);

}
//----------------------------------------------------------------------
