//
//    GRBShock: Class that describes the a Shock between 2 Shells
//    Authors: Nicola Omodei & Johann Cohen Tanugi 
//
#include <iostream>

#include "GRBShock.h"
#include "GRBShell.h"

#define DEBUG 0

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
  M   = sqrt(m1*m1+m2*m2+2.0*m1*m2*(g1*g2-sqrt(g1*g1-1.0)*sqrt(g2*g2-1.0)))-(mf);
  gf  = g1;//(m1*g1 + m2*g2)/(mf+M);
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
  drf = (dr2 + dr1)/2.;//cm
  MS = new GRBShell(gf,rf,drf,ef);
  
  b1  = S1->GetBeta();
  b2  = S2->GetBeta();
  b12 = (b1-b2)/(1.0-b1*b2);
  
  double g12 = 1.0/sqrt(1.0 - b12*b12);
  
  double n1  = S1->GetComPartDens(); //1/cm^3
  double n2  = S2->GetComPartDens(); //1/cm^3

  double x = n1/n2;
  gsh = TMath::Max(1.0,pow(x,1./4.)*sqrt(g12));
  nsh  = 4.0 * n1 * gsh;//1/cm^3
  
  B    = 0.39 * sqrt(ab * nsh);
  ts0  = (5.1e8)/pow(B,2.0)/(2.0*gf); //oss
  //   tcom = 2 gf toss 	   
  // Cooling & time scales in observer frame:
  // a) Shock crossing time scale:
  tc    =   rf/(2.* c * gf * gf);  // oss;//drf/(2.* c * gf);//rf/(c*gf)/(2.*gf);// Dafne & Jonathan //oss
  // b) angular time scale
  ta    =   rf/(2.* c * gf * gf);  // oss
  // c) cooling time scale:
  //characteistic gamma factors:
  // a) Minimum:
  gem  = TMath::Max(1.0, (p-2.0)/(p-1.0)* mpc2/mec2 * ae);
  // b) Cooling :
  gec  = TMath::Max(1.0,ts0/tc);//mec2/(Psyn0 * tc);*tc
  // c) Maximum:
  geM  = TMath::Max(1.,1.17e8/sqrt(B)); 
  
  if(DEBUG)
    {
      if(gem<gec)
	std::cout<<" SLOW COOLING (gem = "<<gem<<" gec = "<<gec<<")"<<std::endl;
      else 
	std::cout<<" FAST COOLING (gem = "<<gem<<" gec = "<<gec<<")"<<std::endl;
    }
}

void GRBShock::SetTime(double time)
{
  tsh=time;
}

Double_t GRBShock::Peak(Double_t time,Double_t energy)
{
  Double_t to  = time-tsh;   // oss
  double ts    = ts0 / GammaElectron(energy); //oss
  double width = sqrt(ta*ta + tc*tc + ts*ts);
  return  eff*exp(-pow(fabs(to)/width,1.0));
  
  /*
    
  
  Double_t to  = time-tsh;   // oss
  if(to<=0.0) return 0.0;
  //  return 1.0;
  //tc=ta;
  //Double_t tts = ts0/sqrt(energy);
  Double_t tta = ta;//sqrt(ta*ta+tts*tts);//sqrt(energy);
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
  */
}
//////////////////////////////////////////////////
Double_t GRBShock::SynSpectrum(Double_t time, Double_t energy)
{
  //  if(time<=tsh) return 0.0;
  double to = time-tsh;

  double tcom = 2.0*to*gf;

  double ec,em,eM,gm,gM,gc;
  // meszaros: t
  //  gm = gt(gem,0.0);//tcom);//tcom-2.0*gf*tc);//tcom);//-2.0*gf*tc);
  //  gc = gt(gem,fabs(tcom));
  //  gM = gt(geM,0.0);//tcom);//tcom-2.0*gf*tc);//0);//tcom);//-2.0*gf*tc);
  em = EsynObs(gem);
  ec = EsynObs(gec);
  eM = EsynObs(geM);
  
  Double_t Fv;
  /*  if(energy<=em)
    {
      Fv = pow(energy/em, 1./3.);
    }
  else
    {
      Fv = pow(energy/em,-p/2.)*exp(-energy/(eM));
    }
  
  return Fv; //kev/kev/s
  */
  
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
  return Fv*exp(-energy/eM);
  
}


//////////////////////////////////////////////////
Double_t GRBShock::ICSpectrum(Double_t time, Double_t energy)
{ 
  //  if(time<=tsh) return 0.0;
  //  double to = time-tsh;
  return SynSpectrum(time,energy/(gem*gem)) / (gem*gem)*exp(-energy/(gem*gf*mec2*1.0e3));
}

double GRBShock::ComputeFlux(double time, double energy)
{
  //this has to return Nv in [ph/(cm² s keV)]

  double Nv =Peak(time,energy)* (SynSpectrum(time,energy) + m_IC * ICSpectrum(time,energy))/energy;

  return Nv;
}

void GRBShock::Print()
{
  std::cout<<"--------------------------------------------------"<<std::endl;
  std::cout<<"Final Gamma :"<<gf<<" Mass :"<<mf<<" Radius = "<< rf<<" Whidth: "<<drf<<std::endl;
  std::cout<<"t shock = "<<tsh<<" ta = "<<ta<<" tc = "<<tc<<" ts = "<<ts0<<" Duration = "<<GetDuration()<<std::endl;
  std::cout<<"Efficiency: "<<eff<<" Initial E = "<<ei<<" Final E = "<<ef<<std::endl;
  std::cout <<" Int E (obs): "<<eint_o<<" (com): "<<eint_c<<" B Field = "<<B<<std::endl;
  std::cout<<" Nsh = "<<nsh<<" Gsh = "<<gsh<<std::endl;
  std::cout<<"Gmin = "<<gem<<" Gcool = "<<gec<<" Gmax = "<<geM<<std::endl;
  std::cout<<"Es0 = "<<EsynObs(1)<<" Esm = "<<EsynObs(gem)<<" Esc = "<<EsynObs(gec)<<" EsM = "<<EsynObs(geM)<<" IC/Syn = "<<m_IC<<std::endl;

  std::cout<<"--------------------------------------------------"<<std::endl;
}
