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
  GRBShock(GRBShell *SF, GRBShell *SB, double tshock);
  ~GRBShock()
    {
      delete MS;
    }
  
  GRBShell *MergedShell() {return MS;}
  void SetTime(double time);
  inline void SetSpectralParameters(double alpha, double beta)
    {
      a = alpha;
      b = beta;
    }
  inline   void SetICComponent(double ic){m_IC = ic;}
  
  double GetTime(){return tsh;}
  double GetDuration(){return sqrt(ta*ta + tc*tc + ts*ts);}
  double GetEfficiency(){return eff;}
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
      if(gec/tc > gem) return Es0*gec*gec/(tc*tc);
      return Es0*gem*gem;
    }
  
 private:
  GRBShell *MS; 
  double tsh;
  double eff;
  double gf,ei,ef,mf, eint_o,eint_c,rf,drf;
  double nsh,gsh,B;
  double ta,tc,ts;
  double Es0;
  double gem,gec,geM;
  double a,b,p;
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
