/**
 * @file IRB_routines.cxx
 * @brief
 *
 * Models implemented by Luis C Reyes (lreyes04@calpoly.edu) and Jim Chiang (jchiang@slac.stanford.edu)
 * Refer to eblAtten/data/README_EBLmodels.txt for more information
 * 
 * $Header$
 */


//-------------------------------------------------------------------

#include <cmath>
#include <fstream>
#include <iostream>

#include "facilities/commonUtilities.h"

#include "AsciiTableModel.h"
#include "Primack05.h"

using namespace std;

namespace IRB {

float calcKneiske(float energy, float redshift);
float calcPrimack05(float energy, float redshift);
float calcKneiske_HighUV(float energy, float redshift);
float calcStecker05(float energy, float redshift);
float calcStecker05_FE(float energy, float redshift);
float calcFranceschini(float energy, float redshift);
float calcFinke(float energy, float redshift);
float calcSalamonStecker(float energy, float redshift);
float calcGeneric (float energy, float redshift);
float calcGilmore09(float energy, float redshift);
float calcGilmore12_fixed(float energy, float redshift);
float calcGilmore12_fiducial(float energy, float redshift);
float calcInoue13(float energy, float redshift);
float calcDominguez11(float energy, float redshift);
float calcScully14_highOp(float energy, float redshift);
float calcScully14_lowOp(float energy, float redshift);
float calcKneiskeDole10(float energy, float redshift);
float calcKneiskeDole10_CMB(float energy, float redshift);
float calcHelgasonKashlinsky12(float energy, float redshift);

float calcHelgasonKashlinsky12(float energy, float redshift) {
   static std::string datadir(facilities::commonUtilities::getDataPath("eblAtten"));
   static std::string od_file =
      facilities::commonUtilities::joinPath(datadir, "opdep_HELGASONandKASHLINSKY_2012.dat");
   
   static AsciiTableModel HelgasonKashlinsky12(od_file);
   return HelgasonKashlinsky12.value(energy, redshift);
}


float calcKneiskeDole10_CMB(float energy, float redshift) {
   static std::string datadir(facilities::commonUtilities::getDataPath("eblAtten"));
   static std::string od_file =
      facilities::commonUtilities::joinPath(datadir, "opdep_KNEISKEandDOLE_2010.dat");

   static AsciiTableModel kneiskedole10_cmb(od_file); 
   return kneiskedole10_cmb.value(energy, redshift);
}


float calcKneiskeDole10(float energy, float redshift) {
   static std::string datadir(facilities::commonUtilities::getDataPath("eblAtten"));
   static std::string od_file =
      facilities::commonUtilities::joinPath(datadir, "opdep_KNEISKEandDOLE_2010_noCMB.dat");
   
   static AsciiTableModel kneiskedole10(od_file);
   return kneiskedole10.value(energy, redshift);
}

float calcScully14_highOp(float energy, float redshift) {
   static std::string datadir(facilities::commonUtilities::getDataPath("eblAtten"));
   static std::string od_file =
      facilities::commonUtilities::joinPath(datadir, "opdep_SCULLYetal2014_highOp.dat");
   
   static AsciiTableModel scully14_highOp(od_file);
   return scully14_highOp.value(energy, redshift);
}

float calcScully14_lowOp(float energy, float redshift) {
   static std::string datadir(facilities::commonUtilities::getDataPath("eblAtten"));
   static std::string od_file =
      facilities::commonUtilities::joinPath(datadir, "opdep_SCULLYetal2014_lowOp.dat");

   static AsciiTableModel scully14_lowOp(od_file);
   return scully14_lowOp.value(energy, redshift);
}


float calcDominguez11(float energy, float redshift) {
   static std::string datadir(facilities::commonUtilities::getDataPath("eblAtten"));
   static std::string od_file =
      facilities::commonUtilities::joinPath(datadir, "opdep_DOMINGUEZetal_2011.dat");
   
   static AsciiTableModel dominguez11(od_file);
   return dominguez11.value(energy, redshift);
}

float calcInoue13(float energy, float redshift) {
   static std::string datadir(facilities::commonUtilities::getDataPath("eblAtten"));
   static std::string od_file =
      facilities::commonUtilities::joinPath(datadir, "opdep_INOUEetal_2013.dat");
   
   static AsciiTableModel inoue13(od_file);
   return inoue13.value(energy, redshift);
}

float calcGilmore09(float energy, float redshift) {
   static std::string datadir(facilities::commonUtilities::getDataPath("eblAtten"));
   static std::string od_file =
      facilities::commonUtilities::joinPath(datadir, "opdep_GILMOREetal_2009.dat");
   
   static AsciiTableModel gilmore09(od_file);
   return gilmore09.value(energy, redshift);
}

float calcGilmore12_fixed(float energy, float redshift) {
   static std::string datadir(facilities::commonUtilities::getDataPath("eblAtten"));
   static std::string od_file = 
      facilities::commonUtilities::joinPath(datadir, "opdep_fixed_Gilmore2012.dat");

   static AsciiTableModel gilmore_fixed(od_file);
   /// The Gilmore et al 2012 tables use MeV, so convert back from GeV.
   return gilmore_fixed.value(energy*1e3, redshift);
}

float calcGilmore12_fiducial(float energy, float redshift) {
   static std::string datadir(facilities::commonUtilities::getDataPath("eblAtten"));
   static std::string od_file = 
      facilities::commonUtilities::joinPath(datadir, "opdep_fiducial_Gilmore2012.dat");
   static AsciiTableModel gilmore_fiducial(od_file);
   /// The Gilmore et al 2012 tables use MeV, so convert back from GeV.
   return gilmore_fiducial.value(energy*1e3, redshift);
}


float calcSalamonStecker(float energy, float redshift) {
   static std::string datadir(facilities::commonUtilities::getDataPath("eblAtten"));
   static std::string od_file =
      facilities::commonUtilities::joinPath(datadir, "opdep_SalamonStecker98.dat");
   static AsciiTableModel salamon_stecker(od_file);
   return salamon_stecker.value(energy, redshift);
}

float calcFinke(float energy, float redshift) {
   static std::string datadir(facilities::commonUtilities::getDataPath("eblAtten"));
   static std::string od_file =
      facilities::commonUtilities::joinPath(datadir, "opdep_Finke_2009.dat");
   static AsciiTableModel finke(od_file);
   return finke.value(energy, redshift);  
}

float calcFranceschini(float energy, float redshift) {
   static std::string datadir(facilities::commonUtilities::getDataPath("eblAtten"));
   static std::string od_file =
      facilities::commonUtilities::joinPath(datadir, "opdep_Franceschini_2008.dat");
   static AsciiTableModel franceschini(od_file);
   return franceschini.value(energy, redshift);
}

float calcKneiske(float energy, float redshift) {
   static std::string datadir(facilities::commonUtilities::getDataPath("eblAtten"));
   static std::string od_file =
      facilities::commonUtilities::joinPath(datadir, "opdep_KNEISKEetal2004_bestfit.dat");
   static AsciiTableModel kneiske(od_file);
   return kneiske.value(energy, redshift);
}


float calcKneiske_HighUV(float energy, float redshift) {
   static std::string datadir(facilities::commonUtilities::getDataPath("eblAtten"));
   static std::string od_file =
      facilities::commonUtilities::joinPath(datadir, "opdep_KNEISKEetal2004_highUV.dat");
   static AsciiTableModel kneiske_uv(od_file);
   return kneiske_uv.value(energy, redshift);
}


float calcPrimack05(float energy, float redshift) {
   return Primack05::instance().value(energy, redshift);
}


float calcGeneric (float energy, float redshift){

//Parametric representation of tau(E,z) from Justin Finke


double E_0 = 80./(1.+redshift);    //energy in units of GeV
double E_1 = 450./(1.+redshift);
double E = (double)energy;

float tau;


if (E < 10. || redshift <= 0.) return 0.;

if (E < E_0)
   tau = (float)(2.6 * (1.+redshift) * pow(E_0/E_1, 2.) * pow(E/E_0, 4.));
   else if (energy < E_1)
        tau = (float)(2.6 * (1.+redshift) * pow(E/E_1, 2.));
		else
		   tau = (float)(2.6 * (1+redshift) * pow(E/E_1, 0.6));
		   
		 
		   
 return tau;
 
 }
		   


float calcStecker05(float energy, float redshift){
//EBL model 3:  Stecker, Malkan, and Scully
//              Astro-ph 0510449
//		Valid for opacities  0.01 < tau < 100


double tau1, tau2, tauvalue;

int zindex=0, MAXZINDEX=9;

double zvalue[9] = {0., 0.03, 0.117, 0.2, 0.5, 1., 2., 3., 5.};

double EMIN [9] = {80., 60., 35., 25., 15., 10., 7., 5., 4.};


double coeff[9][5] = {{0., 0., 0., 0., 0.},
                     {-0.020228, 1.28458, -29.1498, 285.131, -1024.64},
                     {0.010677, -0.238895, -1.004, 54.1465, -313.486},
		     {0.0251369, -0.932664, 11.4876, -45.9286, -12.1116},
		     {-0.0221285, 1.31079, -28.2156, 264.368, -914.546},
		     {-0.175348, 8.42014, -151.421, 1209.13, -3617.51},
		     {-0.311617, 14.5034, -252.81, 1956.45, -5671.36},
		     {-0.34995, 16.0968, -277.315, 2121.16, -6077.41},
		     {-0.321182, 14.6436, -250.109, 1897.00, -5390.55}};



if (redshift < 0.){
   std::cerr<<"Invalid redshift (z<0) ..."<<std::endl;
   redshift=0.;
   } else if (redshift > 5.){
      std::cerr<<"This model is only valid for z <= 5."<<std::endl;
      redshift=5.;
      }

double x = log10(energy*1e+09);


//Find zindex
 for(int i=0; i<MAXZINDEX-1; i++)
    if(redshift >= zvalue[i] && redshift < zvalue[i+1]) zindex = i;
  if(redshift >= zvalue[MAXZINDEX-1]) zindex = MAXZINDEX-1;


if (energy <= EMIN[zindex])
   return 0.;

if (redshift == 0.0)
   return 0.;

if (zindex == 0){
   tau1 = 0.;
   tau2 = coeff[1][0]*pow(x,4.)+coeff[1][1]*pow(x, 3.)+coeff[1][2]*x*x+coeff[1][3]*x+coeff[1][4];
   tauvalue =tau2*redshift/zvalue[1];
   }
   else if (zindex < MAXZINDEX-1){
   tau1 = coeff[zindex][0]*pow(x,4.)+coeff[zindex][1]*pow(x, 3.)
                        +coeff[zindex][2]*x*x+coeff[zindex][3]*x+coeff[zindex][4];
   tau2 = coeff[zindex+1][0]*pow(x,4.)+coeff[zindex+1][1]*pow(x, 3.)
                        +coeff[zindex+1][2]*x*x+coeff[zindex+1][3]*x+coeff[zindex+1][4];
   tauvalue= tau1 + (tau2-tau1)*(redshift-zvalue[zindex])/(zvalue[zindex+1]-zvalue[zindex]);

   } else
      tauvalue = coeff[MAXZINDEX-1][0]*pow(x,4.)+coeff[MAXZINDEX-1][1]*pow(x, 3.)
                        +coeff[MAXZINDEX-1][2]*x*x+coeff[MAXZINDEX-1][3]*x+coeff[MAXZINDEX-1][4];


return pow(10., tauvalue);

}


float calcStecker05_FE(float energy, float redshift){
//EBL model 3:  Stecker, Malkan, and Scully
//              Astro-ph 0510449
//		Valid for opacities  0.01 < tau < 100


double tau1, tau2, tauvalue;

int zindex=0, MAXZINDEX=9;

double zvalue[9] = {0., 0.03, 0.117, 0.2, 0.5, 1., 2., 3., 5.};

double EMIN [9] = {80., 60., 35., 25., 15., 10., 7., 5., 4.};

double coeff[9][5] = {{0., 0., 0., 0., 0.},
		     {-0.020753, 1.31035, -29.6157, 288.807, -1035.21},
			 {0.022352, -0.796354, 8.95845, -24.8304, -79.0409},
		     {0.0258699, -0.960562, 11.8614, -47.9214, -8.90869},
		     {0.0241367, -0.912879, 11.7893, -54.9018, 39.2521},
		     {-0.210116, 10.0006, -178.308, 1412.01, -4190.38},
		     {-0.397521, 18.3389, -316.916, 2431.84, -6991.04},
		     {-0.344304, 15.8698, -273.942, 2099.29, -6025.38},
		     {-0.28918, 13.2673, -227.968, 1739.11, -4969.32}};


if (redshift < 0.){
   std::cerr<<"Invalid redshift (z<0) ..."<<std::endl;
   redshift=0.;
   } else if (redshift > 5.){
      std::cerr<<"This model is only valid for z <= 5."<<std::endl;
      redshift=5.;
      }

double x = log10(energy*1e+09);


//Find zindex
 for(int i=0; i<MAXZINDEX-1; i++)
    if(redshift >= zvalue[i] && redshift < zvalue[i+1]) zindex = i;
  if(redshift >= zvalue[MAXZINDEX-1]) zindex = MAXZINDEX-1;


if (energy <= EMIN[zindex])
   return 0.;

if (redshift == 0.0)
   return 0.;

if (zindex == 0){
   tau1 = 0.;
   tau2 = coeff[1][0]*pow(x,4.)+coeff[1][1]*pow(x, 3.)+coeff[1][2]*x*x+coeff[1][3]*x+coeff[1][4];
   tauvalue =tau2*redshift/zvalue[1];
   }
   else if (zindex < MAXZINDEX-1){
   tau1 = coeff[zindex][0]*pow(x,4.)+coeff[zindex][1]*pow(x, 3.)
                        +coeff[zindex][2]*x*x+coeff[zindex][3]*x+coeff[zindex][4];
   tau2 = coeff[zindex+1][0]*pow(x,4.)+coeff[zindex+1][1]*pow(x, 3.)
                        +coeff[zindex+1][2]*x*x+coeff[zindex+1][3]*x+coeff[zindex+1][4];
   tauvalue= tau1 + (tau2-tau1)*(redshift-zvalue[zindex])/(zvalue[zindex+1]-zvalue[zindex]);

   } else
      tauvalue = coeff[MAXZINDEX-1][0]*pow(x,4.)+coeff[MAXZINDEX-1][1]*pow(x, 3.)
                        +coeff[MAXZINDEX-1][2]*x*x+coeff[MAXZINDEX-1][3]*x+coeff[MAXZINDEX-1][4];


return pow(10., tauvalue);

}

} // namespace IRB
