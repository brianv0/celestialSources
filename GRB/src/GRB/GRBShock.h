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
  ~GRBShock(){;}
  GRBShell *MergedShell() {return MS;}
  void SetTime(double time);
  double GetTime(){return tsh;}
  double Peak(double time, double energy);
  double SynSpectrum(double energy);
  double ICSpectrum(double energy);
  double ComputeFlux(double time, double energy);
  
 private:
  GRBShell *MS; 
  double tsh;
  double gf,ei,ef,mf, eint_o,eint_c,rf,drf;
  double ta,tc;
  double Esyn,Eic;
  float a,b;
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
