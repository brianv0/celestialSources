//
//    GRBShock: Class that describes the a Shock between 2 Shells
//    Authors: Nicola Omodei & Johann Cohen Tanugi 
//
#include <iostream>

#include "GRBShock.h"
#include "GRBShell.h"

using namespace cst;
//////////////////////////////////////////////////
// S2 --> S1 ->
GRBShock::GRBShock(GRBShell *S1, GRBShell *S2, double tshock)
{
  p = 2.5;
  m_IC =  1.0;
  
  tsh = tshock;
  
  double g1 = S1->GetGamma();
  double g2 = S2->GetGamma();
  
  double r1 = S1->GetRadius();
  double r2 = S2->GetRadius();

  double dr1 = S1->GetThickness();
  double dr2 = S2->GetThickness();

  double m1 = S1->GetMass();
  double m2 = S2->GetMass();
  
  double e1 = S1->GetEnergy();
  double e2 = S2->GetEnergy();
  
  // Initial energy of the final shell: (obs)
  ei  = (m1*g1+m2*g2) * c2; //erg

  // mass of the final shell:
  mf  = m1+m2; //g
  //////////////////////////////////////////////////  
  // my model:
  double M,b12,b1,b2,nc,Psyn;
  
  //  std::cout<<"--------------------My Model--------------------"<<std::endl;
  M   = sqrt(m1*m1+m2*m2+2.0*m1*m2*(g1*g2-sqrt(g1*g1-1.0)*sqrt(g2*g2-1.0)))-(mf);//g
  
  gf  = (m1*g1+m2*g2)/(mf+M);
  // Final: (obs) 
  ef  = mf * gf * c2;       //erg
  // Internal : (obs) 
  eint_o = ei - ef;         //erg
  eint_c = (eint_o)/gf;     //erg
  eff = eint_o/ei;
  //////////////////////////////////////////////////
  // Radius of the shell:
  rf  = r2;//cm

  // Thickness of the final shell
  drf = (dr2 + dr1);///2.;//cm

  MS = new GRBShell(gf,rf,drf,ef,mf);
  
  b1  = S1->GetBeta();
  b2  = S2->GetBeta();
  b12 = (b1-b2)/(1-b1*b2);
  nc  = S1->GetComPartDens(); //1/cm^3
  gsh = 1.0/sqrt(1.0 - b12*b12);
  nsh = 4.0 * nc * gsh;//1/cm^3
  B   = 0.39 * sqrt(ab * nsh) * gsh;//G

  double Psyn0 = 6.6e-7 * pow(B,2.); //KeV/s (/g^2) 
  Es0 =1.1447e-11 * B * gf; //keV (/g^2)

  //  std::cout<<mec2/Psyn0<<std::endl;

  // Cooling & time scales in observer frame:
  // a) Hydrodinamical time scale:
  tc    =   drf/c;//(2.* gf * c);  // oss
  // b) angular time scale
  ta    =   rf/(2.* c * gf * gf);  // oss
  // c) cooling time scale:
  ts    = (mec2/Psyn0)/(2.*gf)*sqrt(Es0); //oss //*sqrt(Es0); // s (ts(E) ~ 1/sqrt(E[kev]) 
  //characteristic gamma factors:
  // a) Minimum:
  gem  =   TMath::Max(1.0, mpc2/mec2 * gsh * ae);
  // b) Cooling :
  gec  = TMath::Max(1.0,ts*(2.*gf)*sqrt(Es0));//mec2/(Psyn0 * tc);*tc
  // c) Maximum:
  geM  = TMath::Max(1.,1.17e8/sqrt(B)); 
  
  // characteristic synchrotron energies
  /*
    Esm = Es0 * gem * gem;  // keV
    Esc = Es0;// * gec * gec;  // keV 
    EsM = Es0 * geM * geM;  // keV
  */  
  // IC
  //  Eicm = gem*gem * Esm; 
  // EicM = geM*geM * EsM;

  
  /*
    std::cout<<"gf = "<<gf<<" Esyn m, M= "<<Esm<<" "<<EsM<<" Eic m, M = "<<Eicm<<" "<<EicM<<" ts(100kev) = "<<ts <<std::endl;
    std::cout<<"tc,ta,ts = "<<tc<<" "<<ta<<" "<<ts<<std::endl;
  */
}

void GRBShock::SetTime(double time)
{
  tsh=time;
}

Double_t GRBShock::Peak(Double_t time,Double_t energy)
{
  
  Double_t to  = time-tsh;   // oss
  if(to<=0.0) return 0.0;
  //  return 1.0;

  Double_t tts = ts;
  Double_t tta = sqrt(ta*ta+tts*tts)/sqrt(energy);
  //Double_t tta = ta;
  
  Double_t F1 = 1./pow(1.+       to/tta,2.);
  Double_t F2 = 1./pow(1.+(to - tc)/tta,2.);
  const Double_t H = eint_o /1.0e52;

  if(to <= tc)
    {
      return H*(1.- F1);
    }
  else
    {
      return H*(F2 - F1);
    }
}
//////////////////////////////////////////////////
Double_t GRBShock::SynSpectrum(Double_t time, Double_t energy)
{
  if(time<=tsh) return 0.0;
  double tt = time-tsh;
  //  return 1.0;
  double ec,em,gc;
  /*
  if(tt<tc) 
    {
      gc = TMath::Max(1.0,ts/tt); //?    
      gm = TMath::Max(1.0,gem);
    }
  else
    {
      gc = TMath::Max(1.0,ts/(tc+tt));
      gm = TMath::Max(1.0,gem*ts/(ts+tt*gem));
    }
  */
  gc = TMath::Max(1.0,gec/tc); //?    

  //  double ec = Esc*TMath::Max(1.0,pow(ts/(time-tsh),2.0)); //?
  em = Es0*gem*gem;
  ec = Es0*gc*gc;
  Double_t Fv;
  if(ec <= em) //FAST COOLOING REGIME
    {
      if(energy<=ec)
	Fv = pow(energy/ec,1./3.);
      else if(energy<em)
	Fv = pow(energy/ec,-1./2.);
      else
	Fv = pow(em/ec,-1./2.)*pow(energy/em,-p/2.);
    }
  else  // SLOW COOLING REGIME:
    {
      if(energy<=em)
	Fv = pow(energy/em,1./3.);
      else if(energy<ec)
	Fv = pow(energy/em,-(p-1.)/2.);
      else
	Fv = pow(ec/em,-(p-1.)/2.)*pow(energy/ec,-p/2.);
    }
  return Fv*exp(-energy/(Es0*geM*geM));
}


//////////////////////////////////////////////////
Double_t GRBShock::ICSpectrum(Double_t time, Double_t energy)
{ 
  if(time<=tsh) return 0.0;
  double tt = time-tsh;
  /*
  if(tt<tc) 
    {
    gem = TMath::Max(1.0,gem);
    }
  else
  {
  gem = TMath::Max(1.0,gem*ts/(ts+tt*gem));
  }
  */
  return SynSpectrum(time,energy/(gem*gem)) / (gem*gem);
}

/*
  //  Double_t N;
  Double_t E;
  if(energy > EicM) return 0.0;   
  if (Eicm>emin)
    {
      //      N = (pow(Eicm  ,a+1.)- pow(emin,a+1.))/((a+1.)*pow(Eicm,a)) +
      //	(pow(emax,b+1.)- pow(  Eicm,b+1.))/((b+1.)*pow(Eicm,b));
      
      E = (pow(Eicm  ,a+2.)- pow(emin,a+2.))/((a+2.)*pow(Eicm,a)) +
	(pow(emax,b+2.)- pow(  Eicm,b+2.))/((b+2.)*pow(Eicm,b));
    }
  else
    {
      //      N = 
      //	(pow(emax,b+1.)- pow(  Eicm,b+1.))/((b+1.)*pow(Eicm,b));
      E = 
	(pow(emax,b+2.)- pow(  Eicm,b+2.))/((b+2.)*pow(Eicm,b));
    }
  
  if(energy < Eicm) 
    return pow(energy/Eicm,a)/E ;
  return pow(energy/Eicm,b)/E;
}
*/

double GRBShock::ComputeFlux(double time, double energy)
{
  return Peak(time,energy) * (SynSpectrum(time,energy) + m_IC * ICSpectrum(time,energy))/energy;
}

void GRBShock::Print()
{
  std::cout<<"--------------------------------------------------"<<std::endl;
  std::cout<<"Final Gamma :"<<gf<<" Mass :"<<mf<<" Radius = "<< rf<<" Whidth: "<<drf<<std::endl;
  
  std::cout<<"t shock = "<<tsh<<" ta = "<<ta<<" tc = "<<tc<<" ts = "<<ts<<std::endl;
  std::cout<<"Efficiency: "<<eff<<" Initial E = "<<ei<<" Final E = "<<ef<<std::endl;
  std::cout <<" Int E (obs): "<<eint_o<<" (com): "<<eint_c<<" B Field = "<<B<<std::endl;
  std::cout<<"Gmin = "<<gem<<" Gcool = "<<gec/tc<<" Gmax = "<<geM<<std::endl;
  std::cout<<"Es0 = "<<Es0<<" Esm = "<<Es0*gem*gem<<" Esc = "<<Es0*gec*gec/(tc*tc)<<" EsM = "<<Es0*geM*geM<<" IC/Syn = "<<m_IC<<std::endl;

  std::cout<<"--------------------------------------------------"<<std::endl;
}
