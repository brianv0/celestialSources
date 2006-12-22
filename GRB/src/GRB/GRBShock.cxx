//
//    GRBShock: Class that describes the a Shock between 2 Shells
//    Authors: Nicola Omodei & Johann Cohen Tanugi 
//

#include <iostream>
#include <math.h>
#include "GRBShock.h"

using namespace std;

GRBShock::GRBShock(GRBShell* Sh1, GRBShell* Sh2) 
{ 
  m_sum=1.0;
  double m1 = Sh1->Mass();
  double m2 = Sh2->Mass();
  
  double g1 = Sh1->Gamma();
  double g2 = Sh2->Gamma();

  double dr1=Sh1->Thickness();
  double dr2=Sh2->Thickness();
  
  //  double b1=sqrt(1.-1./pow(g1,2));
  //  double b2=sqrt(1.-1./pow(g2,2));
  //  double m11=pow(m1,2);
  //  double m22=pow(m2,2);
  //  double m12=m1*m2;

  double eb=cst::alphab/0.01;
  //Number and density of paerticles in the shocked material:
  m_npart =m1*cst::c2/(cst::mpc2);
  double n1=m1*cst::c2/(cst::mpc2*Sh1->VolCom());

  // relative gamma between the 2 shells:
  m_gr=(m1*g1+m2*g2)/(m1+m2);
  // Gamma of the resulting shell: 
  m_gf=sqrt((m1*g1+m2*g2)/(m1/g1+m2/g2));
  // Internal gamma:
  m_gi=m_gr-m_gf;
  //Gamma Shock
  //  m_gsh =1.0+(m_Eint/((m1+m2)*cst::c2));
  m_gsh=sqrt(m_gi/2.0);
  if (m_gsh<1.0) {m_gsh=1.0;}

  //Calculated by Nich
  //  m_Eint=cst::c2*(sqrt(m11+m22+2*m12*(g1*g2)*(1-b1*b2))-(m1+m2));
  //  m_gf=(m1*g1+m2*g2)/sqrt(m11+m22+2*m12*(g1*g2)*(1-b1*b2));
  
  // See Piran 1999
  // In ternal Energy:
  m_Eint=m1*cst::c2*(g1-m_gf)+m2*cst::c2*(g2-m_gf);
  // Value for the resulting Shell:
  Sh1->setMass(m1+m2);
  Sh1->setGamma(m_gf);
  Sh1->setThickness(dr1+dr2);
  Sh1->setRadius(Sh2->Radius());
  //Values of two shocks data members:
  m_radius= Sh1->Radius();
  m_mass= Sh1->Mass();
  
  //Comoving Volume [cm^3]
  m_VolCom = Sh1->VolCom();


  // Piran
  //  m_npart =m_Eint*cst::erg2MeV/(m_gf*cst::mpc2);
  // Omodei
    // m_npart=(m1+m2)*cst::c2*cst::erg2MeV/(cst::mpc2);

  double E52=m_Eint/1e+52;
  //m_npart/m_VolCom;
  double X=0.0;
  double Y=0.0;

  //B Field
  m_Beq=(.23/5.9)*sqrt(eb)*m_gf*sqrt(n1);  
  
  //Gamma Lorenz e-
  m_gemin =cst::alphae*m_gf*(cst::mpc2)/(cst::mec2)*(cst::p-2.0)/(cst::p-1.0);
  m_gemax = 2.36e+8*1/sqrt(1.0+X+Y)*pow(eb,-1/4)*pow(m_gf,-1/2)*pow(n1,-1/4);
  m_gecoo = 5.24e+4/((1+X+Y)*eb)*pow(m_gf,-1/3)*pow(n1,-2/3)/pow(E52,1/3);
  
  //  m_riset=(dr1+dr2)/(cst::c);
  m_riset=(dr1)/(cst::c);
}

double GRBShock::Esyn(double gammae)
{
  return (4.3e+6)*(m_gf*m_Beq*pow(gammae,2))*cst::hplanck;
}

double GRBShock::Eic(double gammae)
{
  return (4.3e+6)*(m_gf*m_Beq*pow(gammae,4))*cst::hplanck;
}

double GRBShock::OpticalDepht()
{  
  return 1.0;
  // cst::st*cst::csi*n1*min(dr,c*tsyn*((esyn/1.d+5)**(-1/2)))
}

double GRBShock::tsyn(double ee)
{
  return 6.22e+11*pow(m_Beq,(-5./2.))*sqrt((m_gf*cst::hplanck)/ee);
}


double GRBShock::fred(double ee,double tt)
{
  double decayt;
  double riset;
  double Emax=Esyn(m_gemax);
  double Emin=Esyn(m_gecoo);
  //riset=(m_riset);//*(1.0+(Emax-ee)/(Emax-Emin)));
  //  cout<<m_riset<<"  "<<Emin<<"  "<<Emax<<"  "<<m_gf<<endl;
  //riset=(m_riset*(m_gf/(ee)));
  riset=(m_riset*sqrt(1.0e+10/pow(ee,2.0)));
  if (riset<5e-2) riset=5e-2;
  //  cout<<riset<<endl;
  double tp=m_tobs+riset;
  
  if (tsyn(ee)>riset){
    decayt=tsyn(ee);
  }  else {
    decayt=riset;
  } 
  double norma=(riset/2.0+decayt);
  //  double norma=decayt+riset-(riset*exp((m_tobs-tp)/riset));
  
  //  cout<<" "<<decayt<<" "<<riset<<" "<<tp<<" "<<norma<<endl;
  if (norma <=0) 
    {
      return 0.0;
    }  
  else if (tt<=m_tobs)
    {
      return 0.0;
    } else if (tt<=tp)
      {
	//	return (exp(-(tp-tt)/riset))/norma;
	return (tt-m_tobs)/(riset*norma);
      } else 
	{
	  return (exp(-(tt-tp)/decayt))/norma;
	}
}

double GRBShock::Fsyn(double ee,double tt)
{
  if (m_sum<=0) return 0.0;
  double fmax=m_Beq*pow(m_gf,2);
  //erg/s/eV
  double flux;
  if (Esyn(m_gemin)>=Esyn(m_gecoo))
    {
      //      cout<<"Fast Cooling"<<endl;
      if (ee<=Esyn(m_gecoo))
	flux = pow(ee/Esyn(m_gecoo),1./3.);
      else if (ee<=Esyn(m_gemin)) 
	flux = pow(ee/Esyn(m_gecoo),-1./2.);
      else if (ee<=Esyn(m_gemax))
	flux = pow(Esyn(m_gemin)/Esyn(m_gecoo),-1./2.)*pow(ee/Esyn(m_gemin),-cst::p/2.);
      else 
	flux = 0.0;
    } else 
      { 
	//	cout<<"Slow Cooling"<<endl;
	if (ee<=Esyn(m_gemin))   
	  flux = pow(ee/Esyn(m_gemin),1/3);
	else if (ee<=Esyn(m_gecoo)) 
	  flux = pow(ee/Esyn(m_gemin),-(cst::p-1)/2);
	else if (ee<=Esyn(m_gemax)) 
	  flux = pow(Esyn(m_gecoo)/Esyn(m_gemin),-(cst::p-1)/2)*
	    pow(ee/Esyn(m_gecoo),-cst::p/2);
	else 
	  flux = 0.0;
      }
  //  cout<<"  "<<fmax<<" "<<flux<<"  "<<fred(ee,tt)<<"  "<<tt<<endl;
  return fmax*flux*fred(ee,tt);
}

double GRBShock::Fic(double ee,double tt)
{
  if (m_sum<=0) return 0.0;
  double fmax=m_Beq*pow(m_gf,2);
  //erg/s/eV
  double flux;
  if (Eic(m_gemin)>=Eic(m_gecoo))
    {
      //      cout<<"Fast Cooling"<<endl;
      if (ee<=Eic(m_gecoo))
	flux = pow(ee/Eic(m_gecoo),1./3.);
      else if (Eic(m_gemin)<Esyn(m_gemax))
	{
	  if (ee<=Eic(m_gemin)) 
	    flux = pow(ee/Eic(m_gecoo),-1./2.);
	  else if (ee<=Esyn(m_gemax)) //note the suppression at Asyn max approx KN energy !! 
	    flux = pow(Eic(m_gemin)/Eic(m_gecoo),-1./2.)*
	      pow(ee/Eic(m_gemin),-cst::p/2.);
	  else
	    flux = 0.0;
	} 
      else 
	{
	  if (ee<=Esyn(m_gemax)) 
	    flux = pow(ee/Eic(m_gecoo),-1./2.);
	  else
	    flux = 0.0;
	}
      
    } 
  else 
    {
      //      cout<<"IC SLOW Cooling"<<endl;
      
      if (ee<=Eic(m_gemin))
	flux = pow(ee/Eic(m_gemin),1/3);
      else if(Eic(m_gecoo)<Esyn(m_gemax)) 
	{
	  if (ee<=Eic(m_gecoo)) 
	    flux = pow(ee/Eic(m_gemin),-(cst::p-1)/2);
	  else if (ee<=Esyn(m_gemax)) //note the suppression at Esyn max approx KN energy !! 
	    flux = pow(Eic(m_gecoo)/Eic(m_gemin),-(cst::p-1)/2)*
	      pow(ee/Eic(m_gecoo),-cst::p/2);
	  else 
	    flux = 0.0;
	} 
      else 
	{
	  if (ee<=Esyn(m_gemax))
	    flux = pow(ee/Eic(m_gemin),-(cst::p-1)/2);
	}
    }
  return fmax*flux*fred(ee/1.0e+5,tt);
}


void GRBShock::Write()
{
  cout<< "--------------------Shock's parameters --------------" << endl;
  //  cout<< "Gamma 1    = "<< Sh1->Gamma() << " Gamma 2 = " << Sh2->Gamma() << " Mass 1   = " << Sh1->Mass() <<  " Mass 1   = " << Sh2->Mass() << endl;
  //  cout<< "Gamma fin. = "<< m_gf << " Mass f. = "<< (Sh1->Mass())+(Sh2->Mass()) << " G Shock = " << m_gsh <<  endl;
  cout<< "  Eint = "<< m_Eint <<"  E_rad ="<<m_sum<< " Beq = " << m_Beq   <<endl;
  cout<< "  G.E. min " << m_gemin << " G.E.max " << m_gemax << " G.E. Cooling " << m_gecoo << endl;
  cout<< "  Esyn. min " << Esyn(m_gemin) << " Esyn.max " << Esyn(m_gemax) << " Esyn. Cooling " << Esyn(m_gecoo) << endl;
  cout<< "  E IC min " << Eic(m_gemin) << " EIC max " << Eic(m_gemax) << " EIC Cooling " << Eic(m_gecoo) << endl;
  cout<< "  Vol Com = "<< m_VolCom <<" num ele "<<m_npart<< " n el = "<< m_npart/m_VolCom << endl;
  cout<<" ******* " <<m_radius<< "** T obs = "<<m_tobs<<endl; 
  cout<<m_riset<<endl;
}


