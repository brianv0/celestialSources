/*!
 * \class GRBShock
 *
 * \brief This class implements the shock physics.
 * 
 * The input shell is the material already shocked, 
 * that contains an excess of energy due to the inelastic collision.
 * This class calculates the magnetic field and the parameters to determine
 * the distribution of the shocked accelerated electrons.
 * 
 * \author Nicola Omodei       nicola.omodei@pi.infn.it 
 * \author Johann Cohen-Tanugi johann.cohen@pi.infn.it
 *
 */
#ifndef GRBSHOCK_HH
#define GRBSHOCK_HH 1

#include "GRBConstants.h"
#include "GRBShell.h"

//////////////////////////////////////////////////
class GRBShock
{
 public:
  GRBShock(GRBShell *SF, GRBShell *SB, double tshock, double p=2.5);
  ~GRBShock()
    {
      delete MS;
    }
  
  GRBShell *MergedShell() {return MS;}
  void SetTime(double time);
  /*
    inline void SetSpectralParameters(double alpha, double beta)
    {
    a = alpha;
    b = beta;
    }
  */
  inline   void SetICComponent(double ic){m_IC = ic;}
  
  double GetTime(){return tsh;}
  inline double EsynCom(double g){ return 1.1447e-11 * B * pow(g,2.0);} //kev
  inline double EsynObs(double g){ return EsynCom(g)* gf;} //kev
  inline double GammaElectron(double Eobs){return TMath::Max(1.0,sqrt(8.736e10*Eobs/(B * gf)));}
  inline double Psyn(double g)   { return 6.6e-7 * pow(B,2.) *pow(g,2.0);}     //KeV/s
  inline double gt(double g0, double tcom)    { return TMath::Max(1.0,g0/(1.0+TMath::Max(0.0,tcom)/ts0));}
  inline   double GetDuration(){return  sqrt(ta*ta + tc*tc + pow(ts0/GammaElectron(20.0),2.));}
  inline double GetEfficiency(){return eff;}
  
  double Peak(double time, double energy);
  double SynSpectrum(double time, double energy);
  double ICSpectrum(double time, double energy);
  double ComputeFlux(double time, double energy);
  void Print();
  inline double GetEintO(){return eint_o;}
  inline double GetEintC(){return eint_c;}
  inline double GetGammaf(){return gf;}
  inline double GetEshock(){return nsh* gsh* cst::mpc2;}
  inline double GetEsyn(){return Es0;}
  inline double GetPeak()
    {
      return Es0*gem*gem;
    }
  
 private:
  GRBShell *MS; 
  double tsh;
  double eff;
  double gf,ei,ef,mf, eint_o,eint_c,rf,drf;
  double nsh,gsh,B;
  double ta,tc;
  double Es0,ts0;
  double gem,gec,geM;
  double m_p;
  double m_IC;
  
};

//////////////////////////////////////////////////
class ShockCmp
{
public:
  bool operator()(GRBShock *S1,GRBShock *S2)
  {
    return (S1->GetTime() < S2->GetTime());    
  }
};

//////////////////////////////////////////////////
#endif
