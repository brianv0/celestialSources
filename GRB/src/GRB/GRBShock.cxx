//
//    GRBShock: Class that describes the a Shock between 2 Shells
//    Authors: Nicola Omodei & Johann Cohen Tanugi 
//
#include "GRBShock.h"
#include "GRBShell.h"

using namespace cst;
//////////////////////////////////////////////////
// S2 --> S1 ->
GRBShock::GRBShock(GRBShell *S1, GRBShell *S2, double tshock)
{
  tsh = tshock;
  
  double g1 = S1->GetGamma();
  double g2 = S2->GetGamma();

  //  double r1 = S1->GetRadius();
  double r2 = S2->GetRadius();

  double dr1 = S1->GetThickness();
  double dr2 = S2->GetThickness();

  double m1 = S1->GetMass();
  double m2 = S2->GetMass();
  
  double e1 = S1->GetEnergy();
  double e2 = S2->GetEnergy();
  
  // mass of the final shell:
  mf  = m1+m2; //g

  // Final lorentz factor of the merged shell:
  gf  = TMath::Max(1.,sqrt((m1*g1+m2*g2)/((m1/g1)+(m2/g2))));

  // Initial energy of the final shell: (obs)
  ei  = (m1*g1+m2*g2) * c2; //erg

  // Final: (obs) 
  ef  = mf * gf * c2;       //erg

  // Internal : (obs) 
  eint_o = ei - ef;         //erg
  eint_c = (eint_o)/gf;     //erg

  // Radius of the shell:
  rf  = r2;

  // Thickness of the final shell
  drf = (dr2 + dr1);

  MS = new GRBShell(gf,rf,drf,ef,mf);
  
  // Comoving Volumes and Densities:
  //double v1c = S1->GetComovingVolume();
  //double v2c = S2->GetComovingVolume();
  //double vfc = MS->GetComovingVolume();
  
  //  double n1c = S1->GetComPartDens(); //1/cm^3
  // double n2c = S2->GetComPartDens(); //1/cm^3
  //double nfc = MS->GetComPartDens(); //1/cm^3

  //////////////////////////////////////////////////

  //  double d1 = e1/S1->GetVolume();
  // double d2 = e2/S2->GetVolume();

  double G = g2/g1;
  double Y = e2/e1;
  
  //double g = sqrt(G);
  double y = TMath::Max(1./G*1.00000000000001,TMath::Min(0.9999999999999*G,sqrt(Y)));  
  
  double tv= rf/c*(G*G-1.)/(2.*g2*g2);
  double C1 = pow(G,3.)/((G-y)*(G*y-1.));
  //double C2 = (G*y-1.)/(G-y);
  
  //  double gsh  = 0.5 * sqrt(g1*g2*C2);
  // double gsh2 = 0.5 * sqrt(C1);
  // double esh  = sqrt(d1*d2)*C1*y;
  //double B    = 9.5e3*pow(ab/0.1,0.5)*pow(e2/1e53,0.5)*pow(g1/100.0,-3);// G
  
  Esyn = 1.6  *pow(ae/csi,2.)*pow(ab/0.1,0.5)*pow(e2/1e53,0.5)*pow(g1/100.0,-2.)*pow(tv,-3./2.)*C1*G*G/(G-y);
  Eic  = 140.0*pow(ae/csi,4.)*pow(ab/0.1,0.5)*pow(e2/1e53,0.5)*pow(g1/100.0,-2.)*pow(tv,-3./2.)*C1*C1*G*G/(G-y);
  

  /*
    double gr  = g1*g2-sqrt((g1*g1-1.)*(g2*g2-1.));
    double gsh = sqrt((gr*gr+1.)/2.);
    double nsh = n1c* (4.*gsh + 3.); //1/cm^3
    double esh = gsh * nsh * mpc2; //MeV/cm^3
    
    double gfs = gf*sqrt((1.0+2.0*gf/g1)/(2.0+gf/g1));
    double grs = gf*sqrt((1.0+2.0*gf/g2)/(2.0+gf/g2));
    double bfs = sqrt(1.0 - 1.0/gfs);
    
    std::cout<<gfs<<" "<<grs<<" "<<bfs<<std::endl;
  */
  
  a = -2./3.; // +- 0.15
  b = -2.5;  // +- 0.07
  // Qty, LESI, HESI
  //  Ne   a     b
  //  Fv   a+1   b+1
  // vFv   a+2   b+2
  
  //  double Ub  = ab*esh/erg2meV; //erg/cm^3
  //  double Ue  = ae*esh; //MeV/cm^3
  //  double B     = sqrt(8.*pi*Ub); //G   
  //  double ge    = ae * gsh * mpc2/mec2;

  //  double Ub  = ab * eint_c / v1c; //erg/cm^3
  //  double Ue  = ae * eint_c / v1c * erg2meV; // MeV/cm^3
  
  /*
    double B     = sqrt(8.* pi* Ub); //G   
    double ge    = Ue/csi / (nfc * mec2);
  
    Esyn   = 0.05 * (gf/300.0) * (B/1000.)*pow(ge/100.0,2.); // kev
    Eic    = 500.0*Esyn*pow(ge/100.0,2.);                // kev
    
    double tau_th = st  * n1c * dr1 * g1;
    double tau_pp = st  * n1c * dr1 * g1;
  */

  tc    =   drf/(2.*c);  // oss
  ta    =   rf/(2.*c*gf*gf);  // oss
  //ts    =    
  
  //////////////////////////////////////////////////  
  /*
    std::cout<<" e1= "<<e1<<" m1= "<<m1<<" g1= "<<g1<<" r1= "<<r1<<" dr1= "<<dr1<<" n1c = "<<n1c<<std::endl;
    std::cout<<" e2= "<<e2<<" m2= "<<m2<<" g2= "<<g2<<" r2= "<<r2<<" dr2= "<<dr2<<" n2c = "<<n2c<<std::endl;
    std::cout<<" ef= "<<ef<<" mf= "<<mf<<" gf= "<<gf<<" rf= "<<rf<<" drf= "<<drf<<" nfc = "<<nfc<<std::endl;
    std::cout<<" ei= "<<ei<<" eint_o = "<<eint_o<<" eint_c= "<<eint_c<<std::endl;
    //std::cout<<" Ub "<<Ub<<" Ue "<<Ue<<" B "<<B<<" ge "<<ge<<std::endl;
    std::cout<<" B "<<B<<" Esyn "<<Esyn<<" Eic "<<Eic<<std::endl;
    //std::cout<<" Tau Thompson  "<<tau_th<<" Tau Pair Prod."<<tau_pp<<std::endl;
    std::cout<<" tang = L/c (obs) = "<<ta<<"tcross = l/c (obs) = "<<tc<<std::endl;
    std::cout<<"G = "<<Y<<" Y = "<<Y<<" "<<" tv = "<<tv<<" C1 = "<<C1<<" C2 = "<<C2<<std::endl;
  */
  
  //////////////////////////////////////////////////
  /*
    ofstream fout;
    fout.open("pippo",ios::app);
    fout<<Esyn<<" "<<Eic<<" "<<B<<" "<<eint_o/1e52<<" "<<ta<<" "<<tc<<"\n";
  */
  //////////////////////////////////////////////////
  
}

void GRBShock::SetTime(double time)
{
  tsh=time;
}

Double_t GRBShock::Peak(Double_t time,Double_t energy)
{
  
  Double_t to = time-tsh;   // oss
  //  Double_t ts = sqrt(energy/Esyn);
  if(to<=0.0) return 0.0;
  
  Double_t F1 = 1./pow(1.+       to/ta,2.);
  Double_t F2 = 1./pow(1.+(to - tc)/ta,2.);
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
Double_t GRBShock::SynSpectrum(Double_t energy)
{
  // This method computes the spectrum of syncrotron emission
  // Units : ph/(cm,A2(B s keV)
  Double_t N,E;
  if (Esyn>emin)
    {
      N = (pow(Esyn  ,a+1.)- pow(emin,a+1.))/((a+1.)*pow(Esyn,a)) +
	(pow(emax,b+1.)- pow(  Esyn,b+1.))/((b+1.)*pow(Esyn,b));
      
      E = (pow(Esyn  ,a+2.)- pow(emin,a+2.))/((a+2.)*pow(Esyn,a)) +
	(pow(emax,b+2.)- pow(  Esyn,b+2.))/((b+2.)*pow(Esyn,b));
    }
  else
    {
      N = 
	(pow(emax,b+1.)- pow(  Esyn,b+1.))/((b+1.)*pow(Esyn,b));
      E = 
	(pow(emax,b+2.)- pow(  Esyn,b+2.))/((b+2.)*pow(Esyn,b));
    }
  
  if(energy < Esyn) 
    return pow(energy/Esyn,a)/E ;
  return pow(energy/Esyn,b)/E;
  
}

//////////////////////////////////////////////////
Double_t GRBShock::ICSpectrum(Double_t energy)
{
  Double_t N,E;
  if (Eic>emin)
    {
      N = (pow(Eic  ,a+1.)- pow(emin,a+1.))/((a+1.)*pow(Eic,a)) +
	(pow(emax,b+1.)- pow(  Eic,b+1.))/((b+1.)*pow(Eic,b));
      
      E = (pow(Eic  ,a+2.)- pow(emin,a+2.))/((a+2.)*pow(Eic,a)) +
	(pow(emax,b+2.)- pow(  Eic,b+2.))/((b+2.)*pow(Eic,b));
    }
  else
    {
      N = 
	(pow(emax,b+1.)- pow(  Eic,b+1.))/((b+1.)*pow(Eic,b));
      E = 
	(pow(emax,b+2.)- pow(  Eic,b+2.))/((b+2.)*pow(Eic,b));
    }
  
  if(energy < Eic) 
    return pow(energy/Eic,a)/E ;
  return pow(energy/Eic,b)/E;
}


double GRBShock::ComputeFlux(double time, double energy)
{
  return Peak(time,energy) * (SynSpectrum(energy) + 1.0* ICSpectrum(energy));
}

