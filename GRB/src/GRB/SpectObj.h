#ifndef SpectObj_H
#define SpectObj_H
#include "TH2D.h"
#include "TH1D.h"
#include "TROOT.h"
#include <vector>

struct photon
{
  double time;
  double energy;
};

class SpectObj
{
 public:
  SpectObj(const TH2D* In_Nv);
  
  ~SpectObj()
    {
      delete spec;
      delete times;
      delete Nv;
    }
  TH1D *Integral_E(double e1, double e2);
  TH1D *Integral_E(int ei1, int ei2); 
  TH1D *Integral_T(double t1, double t2, double e1, double e2);
  TH1D *Integral_T(int ti1, int ti2, int ei1, int ei2);
  TH1D *Integral_T(double t1, double t2, double en = 0.0);
  TH1D *Integral_T(int ti1, int ti2, int ei = 1);
  double Integral_E(TH1* Sp, double e1=0.0, double e2=0.0);
  double Integral_E(TH1* Sp, int ei1, int ei2);
  double Integral_T(TH1* Lc, double t1=0.0, double t2=0.0);
  double Integral_T(TH1* Lc, int ti1, int ti2);
  TH1D *ComputeProbability(double enph);
  TH1D *N(TH1D *EN);
  photon GetPhoton(double t0, double enph);
  
  double flux(double time, double enph);
  double interval(double time, double enph);
  double energy(double time, double enph);
  void SaveParameters(double tstart, std::pair<double,double> direction);
  //////////////////////////////////////////////////
  void ScaleAtBATSE(double fluence);
  double GetFluence(double BL=0.0, double BH=0.0);
  double GetT90(double BL=0.0, double BH=0.0);
  //////////////////////////////////////////////////
  TH1D* GetSpectrum(double t=0.0);
  TH1D* GetTimes(double t=0.0);

 private:
  TH2D* Nv;
  int ne,nt;
  double emin,emax;
  double tmin,tmax,deltat;
  TH1D *spec,*times;
};
#endif
