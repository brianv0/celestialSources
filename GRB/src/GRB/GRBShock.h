/*!
 * \class GRBShock
 *
 * \brief This class implements the shock physics.
 * 
 * Computes the shock dynamics between two shells which are colliding one against the other.
 
 * This class calculates the magnetic field and the parameters to determine
 * the distribution of the shocked accelerated electrons (assumed power law).
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
  /* Constructor
     \param FS is the forward shell (slower shell)
     \param BS is the backward shell (faster shell)
     \params tshock is the shock time as measured from the central engine
     \param p is the power law index for the distribution of accelerated electrons
  */
  GRBShock(GRBShell *FS, GRBShell *BS, double tshock, double p=2.5);
  /// Destructor: delete the merged shell
  ~GRBShock()
    {
      delete MS;
    }
  /// return the merged shell
  GRBShell *MergedShell() {return MS;}
  
  void SetTime(double time);

  inline   void SetICComponent(double ic){m_IC = ic;}
  
  double GetTime(){return tsh;}
  double GetRadius(){return rf;}
  /*! 
    \brief Compute the characteristic synchrotron energy
    
    \param g: Lorentz factor of the emitting electron
    The characteristic synchrotron energy, in the emitter co-moving frame, is then:
    \f$ E_{syn}(\gamma)=1.73e-11\frac{B}{Gauss}\gamma^2\f$ (kev)
  */
  inline double EsynCom(double g){ return 1.73e-11 * B * pow(g,2.0);} //kev
  /// return the  observed synchrotron energy
  inline double EsynObs(double g){ return EsynCom(g)* gf;} //kev
  /// Return the Lorentz factor of an electron given its characteristic synchrotron energy
  inline double GammaElectron(double Eobs){return TMath::Max(1.0,sqrt(5.8e10*Eobs/(B * gf)));}
  /// Return the total emitted power by an electron with Lorentz factor g (KeV/s)
  inline double Psyn(double g)   { return 1.0e-6 * pow(B,2.) *pow(g,2.0);}     
  /// Return the duration of the produced pulse, assuming geometrical lag (angular spreading time) and electron cooling
  inline   double GetDuration(){return  sqrt(ta*ta + tc*tc + pow(ts0/GammaElectron(20.0),2.));}
  /// Return the conversion efficiency of energy by the shock
  inline double GetEfficiency(){return eff;}
  /*!
    Represents the pulse shape as a function of:
    \param time in the observed frame
    \param energy in the emitter frame
    
  */
  double Peak(double time, double energy);
  /*!
    Represents the synchrotron spectrum as a function of:
    \param time in the observed frame
    \param energy in the emitter frame
  */
  double SynSpectrum(double time, double energy);
  /*!
    Represents the Inverse Compton spectrum as a function of:
    \param time in the observed frame
    \param energy in the emitter frame
  */
  double ICSpectrum(double time, double energy);
  /// Compute the flux \f$ N(e) \f$ in units of (\f$ ph/(cm^2~s~keV)\f$)
  double ComputeFlux(double time, double energy);
  /// Print out method (for debugging)
  void Print();
  /// return the internal observed energy
  inline double GetEintO(){return eint_o;}
  /// return the internal co-moving energy
  inline double GetEintC(){return eint_c;}
  /// return the Lorentz factor of the merged shell
  inline double GetGammaf(){return gf;}
  /// Return the energy density after the shock (the shocked energy density)
  inline double GetEshock(){return nsh* gsh* cst::mpc2;}
  /*! 
    Return the ratio between the synchrotron energy and 
    the square of the Lorentz factor of the emitting electron: \f$ \frac{E_{syn}(\gamma)}{\gamma^2}\f$ (kev)
  */
  inline double GetEsyn(){return Es0;}
  // Return the peak energy in the Fast Cooling Regime
  inline double GetPeak()
    {
      return Es0*gem*gem;
    }
  
 private:
  GRBShell *MS;
  double tsh;
  double eff;
  double gf,ei,ef,mf, eint_o,eint_c,rf,drf;
  double nsh,esh,gsh,B;
  double ta,tar,tc;
  double Es0,ts0;
  double gem,gec,geM;
  double m_p;
  double m_IC;
  
};

//////////////////////////////////////////////////
/*!
  \class ShockCmp
  
  \brief This class allows the sorting of the GRBShock vector
  
  It sort the shock vector in ascending order
  
*/


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
