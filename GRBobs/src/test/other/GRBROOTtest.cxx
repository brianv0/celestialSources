#include <iostream>
#include "TROOT.h"
#include "TApplication.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TTree.h"

#include "SpectObj/SpectObj.h"
#include "GRBobs/GRBobsSim.h"
#include "GRBobs/GRBobsConstants.h"

using namespace ObsCst;
double EMIN, EMAX, TMIN, TMAX, DT;
int    TBIN, EBIN;
//////////////////////////////////////////////////
// options 
bool savePlots   = false;
int frame        = 10;
bool movie       = false;
bool bandFit     = false;
bool powerlawFit = false;
bool scaled      = false;
bool noise       = false;
double noiseL    = 0.0;
bool ExtraComponent = false;
int extension;
int nlat;
double Band(double *var, double *par);
void help();
void ScanParameters(int Ngrb);

#define DEBUG 0 

#define GenerationArea 1.00 

double AeffLAT(double MeV, double Z)
{
  const double p0 = 9.8e3;
  const double p1 = 120.0;
  const double p2 = 2.0;
  double pr = 1.25*(Z-0.2);
  if(MeV<30.0 || pr<=0) return 0.0; 
  return pr*p0/(1.0+exp(-p2*log(MeV/p1)));
}

double AeffNaI(double MeV, double Z)
{
  if(MeV > 0.01 && MeV <=1.0)
    return 126.0 * 1.0;
  else 
    return 0.0;
}

double AeffBGO(double MeV, double Z)
{
  if(MeV > 0.15 && MeV <=30.0)
    return 200.0 * 1.0;
  else 
    return 0.0;
}



double Band(double *var, double *par)
{

  double a  = par[0];
  double b  = par[1];
  double E0 = pow(10.,par[2]);
  double NT = pow(10.,par[3]);

  double E   = var[0];
  double C   = pow((a-b)*E0/100.0,a-b)*exp(b-a);
  double H   = (a-b) * E0;
  if(E <= H) 
    return E * E * NT * pow(E/100.,a) * exp(-E/E0);
  return E * E * C * NT * pow(E/100.,b); // ph cm^(-2) s^(-1) keV^(-1)
}

TH2D *Load(char name[100]="grb_65540.root")
{
  TFile *mod = new TFile(name);
  TH2D *Nv = (TH2D*) mod->Get("Nv"); //ph m^(-2) s^(-1) keV^(-1)
  
  TMIN = Nv->GetXaxis()->GetXmin();
  TMAX = Nv->GetXaxis()->GetXmax();
  EMIN = Nv->GetYaxis()->GetXmin();
  EMAX = Nv->GetYaxis()->GetXmax();
  
  TBIN = Nv->GetXaxis()->GetNbins();   
  EBIN = Nv->GetYaxis()->GetNbins();
  DT   = TMAX/(TBIN-1.);

  std::cout<<" Matrix Loaded "<<std::endl; 
  std::cout<<" Bin t = "<<TBIN<<" TMIN = "<<TMIN<<" TMAX = "<<TMAX<<std::endl; 
  std::cout<<" Bin e = "<<EBIN<<" EMIN = "<<EMIN<<" EMAX = "<<EMAX<<std::endl; 
  
  Nv->SetXTitle("Time [s]");
  Nv->SetYTitle("Energy [keV]");
  Nv->SetZTitle("N_{v} [ph/m^{2}/s/keV]");
  Nv->GetXaxis()->SetTitleOffset(1.5);
  Nv->GetYaxis()->SetTitleOffset(1.5);
  Nv->GetZaxis()->SetTitleOffset(1.2);
  Nv->GetXaxis()->CenterTitle();
  Nv->GetYaxis()->CenterTitle();
  Nv->GetZaxis()->CenterTitle();

  // Nv ph/m^2/s/keV
  /*
  for(int t=0;t<TBIN;t++)
    for(int e=0;e<EBIN;e++)
      {
	double ene = Nv->GetYaxis()->GetBinCenter(e+1);
	Nv->SetBinContent(t+1,e+1,pow(ene,-2.0));  
      }
  */
  return Nv;
}

void PlotGRB(double enph = 0,double z=0,char name[100]="grb_65540.root",TString name2="GRB_")
{
  std::cout<<name<<std::endl;
  gStyle->SetCanvasColor(10);
  
  TCanvas *cNv = new TCanvas("cNv","Nv");
  cNv->SetLogy();
  cNv->SetLogz();
  std::cout<<" Loading "<<name<<std::endl;
  TH2D *Nv = Load(name); // [ph/keV/s/m²]

  //////////////////////////////////////////////////
  SpectObj *sp = new SpectObj(Nv,0,z);              //ph
  sp->SetAreaDetector(GenerationArea); //like observation sim
  double T90  = sp->GetT90();
  double T100 = TMAX-TMIN;
  
  //////////////////////////////////////////////////
  
  // Ne =  [ph/keV/s/m²]
  gDirectory->Delete("Ne");
  gDirectory->Delete("eNe");
  gDirectory->Delete("LAT");
  gDirectory->Delete("BGO");
  gDirectory->Delete("NaI");
  gDirectory->Delete("e2Ne");
  
  Nv->ProjectionY("Ne");
  TH1D *Ne = (TH1D*) gDirectory->Get("Ne");
  Ne->GetXaxis()->SetTitleOffset(1.1);
  Ne->Scale(DT*1e-4);// Ne =  [ph/keV/cm²] (xspec-like plot lcounts)

  Ne->SetYTitle("Ne [ph/keV/cm^{2}]");
  Ne->SetXTitle(" Energy[keV] ");
  // eNe = Ne * de: [ph/m^2]
  TH1D *eNe = (TH1D*) Ne->Clone();
  eNe->SetTitle("eNe");
  eNe->SetName("eNe");
  eNe->SetYTitle("eNe [keV/keV/cm^{2}]");  
  // e2Ne = e * Ne * de: [keV/cm^2]
  //////////////////////////////////////////////////
  TH1D *LAT = (TH1D*) Ne->Clone();
  LAT->SetTitle("LAT");
  LAT->SetName("LAT");
  LAT->SetYTitle("LAT [counts]");  
  TH1D *NaI = (TH1D*) Ne->Clone();
  NaI->SetTitle("NaI");
  NaI->SetName("NaI");
  NaI->SetYTitle("NaI [counts]");  
  TH1D *BGO = (TH1D*) Ne->Clone();
  BGO->SetTitle("BGO");
  BGO->SetName("BGO");
  BGO->SetYTitle("BGO [counts]");  

  TH1D *e2Ne = (TH1D*) Ne->Clone();
  e2Ne->SetTitle("Fluxes");
  e2Ne->SetName("e2Ne");
  
  e2Ne->SetYTitle(" Flux ");
  e2Ne->SetStats(0);

  Ne->SetLineColor(2);
  eNe->SetLineColor(3);
  e2Ne->SetLineColor(4);

  Ne->SetLineWidth(2);
  eNe->SetLineWidth(2);
  e2Ne->SetLineWidth(2);
  
  
  for(int i=0; i < EBIN; i++)
    {
      double ne = Ne->GetBinContent(i+1); //ph/keV/cm^2
      double de = Ne->GetBinWidth(i+1);//GetBinCenter(i+1);
      double e = Ne->GetBinCenter(i+1);
      double Aeff_LAT = AeffLAT(e*1e-3,1.0); //cm^2
      double Aeff_NaI = AeffNaI(e*1e-3,1.0); //cm^2
      double Aeff_BGO = AeffBGO(e*1e-3,1.0); //cm^2
      
      eNe->SetBinContent(i+1,e * ne); //de
      LAT->SetBinContent(i+1,e * ne * Aeff_LAT); //de
      BGO->SetBinContent(i+1,e * ne * Aeff_BGO); //de
      NaI->SetBinContent(i+1,e * ne * Aeff_NaI); //de
      
      e2Ne->SetBinContent(i+1,e*e*ne);
    }
  
  TCanvas *animatedCanvas = new TCanvas("animatedCanvas","animatedCanvas",500,500);
  animatedCanvas->SetLeftMargin(0.15);
  animatedCanvas->SetLogx();
  animatedCanvas->SetLogy();
  int i1=0;
  int i2;
  TH1D *AnimatedSpectrum;
  if(movie)
    {
      int l=1;
      for(int i=0; i < TBIN; i++)
	{
	  i2=i;
	  if(i%frame==0)
	    {
	      TString Name="Plot/frame";
	      AnimatedSpectrum =   Nv->ProjectionY("Animated",i1,i2);
	      AnimatedSpectrum->GetXaxis()->SetTitle("Energy [keV]");
	      AnimatedSpectrum->GetYaxis()->SetTitle("e^{2} N(e) [erg/cm^{2}/s]");
	      AnimatedSpectrum->GetXaxis()->SetTitleOffset(1.1);
	      AnimatedSpectrum->GetYaxis()->SetTitleOffset(1.5);
	      AnimatedSpectrum->SetLineColor(2);
	      AnimatedSpectrum->SetLineWidth(2);
	      for(int j=0; j < EBIN; j++)
		{
		  double ne = AnimatedSpectrum->GetBinContent(j+1);
		  double ee = AnimatedSpectrum->GetBinCenter(j+1);
		  AnimatedSpectrum->SetBinContent(j+1,ee*ee*ne/624151.0*1e-7);
		}  
	      AnimatedSpectrum->SetMaximum(1e-3);
	      AnimatedSpectrum->SetMinimum(1e-10);
	      if(l==1)	      
		AnimatedSpectrum->Draw("l");
	      else
		AnimatedSpectrum->DrawCopy("lsame");
	      animatedCanvas->Update();
	      if(savePlots)
		{
		  Name+=l;
		  Name+=extension;
		  Name+=".gif";
		  animatedCanvas->Print(Name);
		  l++;
		}
	      i1=i;
	    }
	}
    }
  
  
  
  //////////////////////////////////////////////////
  cNv->cd();
  Nv->Draw("surf");
  
  bool ExtractPhotons = true;
  if(enph<=0) ExtractPhotons = false;
  else if (enph<EMIN)  enph = EMIN;
  
  //////////////////////////////////////////////////
  TCanvas *clc = new TCanvas("clc","clc",600,800);
  TCanvas *csp = new TCanvas("csp","csp");
 
  csp->SetLogx();
  csp->SetLogy();
  
  clc->Divide(1,3);
  
  TH1D *Lct_GBM = sp->Integral_E(GBM1,GBM2);  // ph
  
  /*
    TH1D *Lct_BATSE1 = sp->Integral_E(10.,30.);  // ph
    TH1D *Lct_BATSE2 = sp->Integral_E(1.0e2,3.0e2);  // ph
    TH1D *Lct_BATSE3 = sp->Integral_E(1.0e4,3e4);  // ph
    TH1D *Lct_BATSE4 = sp->Integral_E(1.0e6,3e6);  // ph
  */
    
  TH1D *Lct_BATSE1 = sp->Integral_E(BATSE1,BATSE2);  // ph
  TH1D *Lct_BATSE2 = sp->Integral_E(BATSE2,BATSE3);  // ph
  TH1D *Lct_BATSE3 = sp->Integral_E(BATSE3,BATSE4);  // ph
  TH1D *Lct_BATSE4 = sp->Integral_E(BATSE4,BATSE5);  // ph
  TH1D *Lct_SWIFT = sp->Integral_E(SWIFT1,SWIFT2);  // ph
    
  TH1D *Lct_LAT = sp->Integral_E(LAT1,LAT2);  // ph
  TH1D *Lct_EXT = sp->Integral_E(enph,EMAX);  // ph
  
  Lct_GBM->SetStats(0);
  Lct_LAT->SetStats(0);
  
  Lct_GBM->SetNameTitle("flux [ph]","flux [ph]");
  Lct_GBM->SetXTitle("Time (s)");
  
  Lct_LAT->SetNameTitle("flux [ph]","flux [ph]");
  Lct_LAT->SetXTitle("Time (s)");
  Lct_LAT->SetYTitle("photons");
  
  //  Lct_BATSE->SetName("flux GRB [ph]");
  Lct_GBM->SetName("flux GBM [ph]");
  Lct_LAT->SetName("flux LAT [ph]"); 
  Lct_EXT->SetName("flux EXT [ph]"); 
 
  Lct_BATSE1->SetLineWidth(1);
  Lct_BATSE2->SetLineWidth(1);
  Lct_BATSE3->SetLineWidth(1);
  Lct_BATSE4->SetLineWidth(1);
  Lct_SWIFT->SetLineWidth(1);
  Lct_GBM->SetLineWidth(1);
  Lct_LAT->SetLineWidth(2);
  Lct_EXT->SetLineWidth(1);
  Lct_EXT->SetLineStyle(2);
  
  Lct_BATSE1->SetLineColor(2);
  Lct_BATSE2->SetLineColor(3);
  Lct_BATSE3->SetLineColor(4);
  Lct_BATSE4->SetLineColor(6);
  Lct_SWIFT->SetLineColor(7);

  Lct_LAT->SetLineColor(4);
  //  Lct_EXT->SetLineColor(6);
 
  csp->cd();
  e2Ne->SetMinimum(.1);
  e2Ne->Draw("l");
  eNe->Draw("samel");
  Ne->Draw("samel");

  LAT->SetMarkerStyle(2);
  BGO->SetMarkerStyle(3);
  NaI->SetMarkerStyle(4);

  LAT->Draw("samep");
  NaI->Draw("samep");
  BGO->Draw("samep");

  //  gDirectory->Delete("band");
  // Fit with the band function
  if(bandFit)
    {
      TF1 band("grb_f",Band,EMIN,1.0e+4,4); 
      band.SetParNames("a","b","Log10(E0)","Log10(Const)");
      band.SetLineStyle(2);
      
      band.SetParLimits(0,-2.0,2.0);
      band.SetParLimits(1,-4.0,-1.0);
      band.SetParLimits(2,1.5,3.0);
      band.SetParLimits(3,-10.0,10.0);
      
      band.SetParameter(0,-0.0);
      band.SetParameter(1,-2.0);
      band.SetParameter(2,2.5);
      band.SetParameter(3,0.0);
      e2Ne->Fit("grb_f","QR+","lsame");
      //Ne->Fit("grb_f","QR+","lsame");
      
      double a=band.GetParameter(0);
      double b=band.GetParameter(1);
      double a_err=band.GetParError(0);
      double b_err=band.GetParError(1);

      double E0=pow(10.,band.GetParameter(2));
      double E0_err_p=pow(10.,band.GetParameter(2)+band.GetParError(2));
      double E0_err_m=pow(10.,band.GetParameter(2)-band.GetParError(2));
      double Ep=(a+2)*E0;
      
      std::cout<<"-- BAND FUNCTION eNe FIT ---------------------------"<<std::endl;
      std::cout<<" a     = "<<a<<" +- "<<a_err<<std::endl;
      std::cout<<" b     = "<<b<<" +- "<<b_err<<std::endl;
      std::cout<<" E0    = "<<E0<<" + "<<E0_err_p<<" - "<<E0_err_m<<std::endl;
      std::cout<<" Ep    = "<<Ep<<std::endl;
      std::cout<<" Const = "<<pow(10.,band.GetParameter(3))<<" ph/cm^2/kev, "<<pow(10.,band.GetParameter(3))/T100<<" ph/cm^2/kev/s"<<std::endl;
      std::cout<<"--------------------------------------------------"<<std::endl;

    }
  
  if(powerlawFit)
    {
      TF1 pl("PL_f","x^2 * 10^([0])*(x)^(-[1])",1.0e+3,1.0e+6);
      
      pl.SetLineStyle(2);
      pl.SetLineColor(2);
      //band.Draw("lsame");
      pl.SetParameters(4.0,-2.25);
      pl.SetParLimits(0,-10,10);
      pl.SetParLimits(1,0,4);
      
      //////////////////////////////////////////////////
      e2Ne->Fit("PL_f","QR+","lsame");
      std::cout<<"---PL FUNTION Nv FIT------------------------------"<<std::endl;
      std::cout<<" K       = "<<pow(10.,pl.GetParameter(0))<<" ph/cm^2/kev, at 1keV "<<pow(10.,pl.GetParameter(0))<<" ph/cm^2/kev/s, at 1keV "<<std::endl;
      std::cout<<" F(  1 MeV) = "<<pl.Eval(1e3)<<" ph/m^2/kev"<<std::endl;
      std::cout<<" F( 10 MeV) = "<<pl.Eval(1e4)<<" ph/m^2/kev"<<std::endl;
      std::cout<<" F(100 MeV) = "<<pl.Eval(1e5)<<" ph/m^2/kev"<<std::endl;
      std::cout<<" F(  1 GeV) = "<<pl.Eval(1e6)<<" ph/m^2/kev"<<std::endl;
      std::cout<<" Index = "<<pl.GetParameter(1)<<std::endl;
      std::cout<<"--------------------------------------------------"<<std::endl;
      /*
      eNe->Fit("PL_f","QR+","lsame");
      std::cout<<"---PL FUNTION eNe FIT------------------------------"<<std::endl;
      std::cout<<" No    = "<<pow(10.,pl.GetParameter(0))<<std::endl;
      std::cout<<" Index = "<<pl.GetParameter(1)<<std::endl;
      std::cout<<"--------------------------------------------------"<<std::endl;
      
      e2Ne->Fit("PL_f","QR+","lsame");
      std::cout<<"---PL FUNTION e2Ne FIT------------------------------"<<std::endl;
      std::cout<<" No    = "<<pow(10.,pl.GetParameter(0))<<std::endl;
      std::cout<<" Index = "<<pl.GetParameter(1)<<std::endl;
      std::cout<<"--------------------------------------------------"<<std::endl;
      */
    }
  TLegend *leg = new TLegend(0.11,0.12,0.37,0.25);
  leg->SetFillStyle(0);
  leg->AddEntry(Ne," N(e)  [ph/keV/cm^{2}] ","lp");
  leg->AddEntry(eNe," e N(e) [keV*(ph/keV/cm^{2})] ","lp");
  leg->AddEntry(e2Ne," e^{2} N(e) [keV^{2}*(ph/keV/cm^{2})] ","lp");
  leg->Draw();
  
  TH1D *Counts   = sp->CloneSpectrum(); //ph
  TH1D *Lc       = sp->CloneTimes();  
  TString title = "Photons over ";
  title+=GenerationArea;
  title+=" m^{2}";
  Lc->SetTitle(title);
  Lc->SetXTitle("Time (s)");
  Lc->SetYTitle("photons");

  if(ExtractPhotons)
    {  
      double time = 0.0;
      double Interval;
      double energy;
      double flux;
      int i = 0; 
      
      Interval = sp->interval(time,enph); // s
      time+=Interval;
      
      while(time < TMAX && time>=TMIN)
	{
	  energy   = sp->energy(time,enph);   // keV
	  flux     = sp->flux(time,enph);     // ph/s/m2
	  Counts->Fill(energy); // ph
	  Lc->Fill(time);       // ph
	  Interval = sp->interval(time,enph); // s
      	  std::cout<<" Time (s)  = "<<time
		   <<" Flux (ph/s/m^2)  = "<<flux
		   <<" energy (MeV) = "<<energy*1e-3
		   <<" Interval (s) = "<<Interval<<std::endl;
	  time+=Interval;	  
	  i++;
	} 
    }
  
  if(Counts->GetEntries() == 0 && ExtractPhotons) 
    {
      std::cout<<" sorry, no extracted photons for t <= "<<TMAX<<std::endl;
    }
  
  csp->cd();
  
  if(ExtractPhotons)
    {
      Counts->SetMaximum(1.5*Ne->GetMaximum());
      Counts->SetStats(1);
      for(int ei=1;ei<=EBIN;ei++)
	{
	  double en = Counts->GetBinCenter(ei);
	  double de = Counts->GetBinWidth(ei);
	  double ph = Counts->GetBinContent(ei);
	  double err = sqrt(ph);
	  Counts->SetBinContent(ei,ph/de/GenerationArea*1e-4); //ph/keV/cm^2
	  Counts->SetBinError(ei,err/de/GenerationArea*1e-4); //ph/keV/cm^2
	}
      Counts->Draw("E1same");
      Counts->SetStats(1);
    }
  //////////////////////////////////////////////////      

  Lct_GBM->Scale(1.0/GenerationArea*1e-4);
  Lct_BATSE1->Scale(1.0/GenerationArea*1e-4);
  Lct_BATSE2->Scale(1.0/GenerationArea*1e-4);
  Lct_BATSE3->Scale(1.0/GenerationArea*1e-4);
  Lct_BATSE4->Scale(1.0/GenerationArea*1e-4);
  Lct_SWIFT->Scale(1.0/GenerationArea*1e-4);

  Lct_GBM->SetYTitle("N_{ph}/cm^{2}"); 
  Lct_SWIFT->SetYTitle("N_{ph}/cm^{2}"); 

  if(scaled)
    {
      double MaxGBM = Lct_GBM->GetMaximum();
      double Max1   = Lct_BATSE1->GetMaximum();
      double Max2   = Lct_BATSE2->GetMaximum();
      double Max3   = Lct_BATSE3->GetMaximum();
      double Max4   = Lct_BATSE4->GetMaximum();
      double Max5   = Lct_SWIFT->GetMaximum();
      
      Lct_GBM->Scale(1./MaxGBM);
      Lct_BATSE1->Scale(1./Max1);
      Lct_BATSE2->Scale(1./Max2);
      Lct_BATSE3->Scale(1./Max3);
      Lct_BATSE4->Scale(1./Max4);
      Lct_SWIFT->Scale(1./Max5);
      Lct_SWIFT->SetYTitle("Normalized counts");  
      Lct_GBM->SetYTitle("Normalized counts");  
    }
  else if (noise)
    {
      TRandom rnd;
      double AGBM     = 125.0*6.0;
      double ABATSE   = 2025.0;
      double ASWIFT_BAT = 2000.0;
      double BGBM     = 200.0;
      double BBATSE1   = 220.0;
      double BBATSE2   = 150.0;
      double BBATSE3   = 120.0;
      double BBATSE4   = 80.0;
      double BSWIFT_BAT = 10000.0;
      
      Lct_GBM->Scale(AGBM);
      Lct_BATSE1->Scale(ABATSE);
      Lct_BATSE2->Scale(ABATSE);
      Lct_BATSE3->Scale(ABATSE);
      Lct_BATSE4->Scale(ABATSE);
      Lct_SWIFT->Scale(ASWIFT_BAT);

      Lct_SWIFT->SetYTitle("Counts");  
      Lct_GBM->SetYTitle("Counts");  
      
      for (int i = 1; i<= TBIN; i++)
	{
	  Lct_GBM->SetBinContent(i,rnd.Poisson(BGBM+Lct_GBM->GetBinContent(i)));
	  Lct_BATSE1->SetBinContent(i,rnd.Poisson(BBATSE1+Lct_BATSE1->GetBinContent(i)));
	  Lct_BATSE2->SetBinContent(i,rnd.Poisson(BBATSE2+Lct_BATSE2->GetBinContent(i)));
	  Lct_BATSE3->SetBinContent(i,rnd.Poisson(BBATSE3+Lct_BATSE3->GetBinContent(i)));
	  Lct_BATSE4->SetBinContent(i,rnd.Poisson(BBATSE4+Lct_BATSE4->GetBinContent(i)));
	  Lct_SWIFT->SetBinContent(i,rnd.Poisson(noiseL * BSWIFT_BAT * DT + Lct_SWIFT->GetBinContent(i)));
	}
    }
  Lct_GBM->SetMinimum(-Lct_GBM->GetMaximum()/10.);
  Lct_LAT->SetMinimum(-Lct_LAT->GetMaximum()/10.);
  if(noise)
    {
      clc->cd(1);
      //      gPad->SetLogx();
      gPad->SetLogy();

      Lct_SWIFT->Draw();
      clc->cd(2);
      Lct_GBM->Draw();
      Lct_BATSE1->Draw("same");
      Lct_BATSE2->Draw("same");
      Lct_BATSE3->Draw("same");
      Lct_BATSE4->Draw("same");
    }
  else
    {
      clc->cd(1);
      Lct_SWIFT->Draw("l");
      clc->cd(2);
      Lct_GBM->Draw("l");
      Lct_BATSE1->Draw("lsame");
      Lct_BATSE2->Draw("lsame");
      Lct_BATSE3->Draw("lsame");
      Lct_BATSE4->Draw("lsame");
    }
  TLegend *legL = new TLegend(0.7,0.7,0.99,0.99);
  TLegend *legL1 = new TLegend(0.7,0.7,0.99,0.99);
  char legname[100];
  sprintf(legname, "BATSE %.0f-%.0f keV",BATSE1,BATSE2);
  legL->AddEntry(Lct_BATSE1,legname);
  sprintf(legname, "BATSE %.0f-%.0f keV",BATSE2,BATSE3);
  legL->AddEntry(Lct_BATSE2,legname);
  sprintf(legname, "BATSE %.0f-%.0f keV",BATSE3,BATSE4);
  legL->AddEntry(Lct_BATSE3,legname);
  sprintf(legname, "BATSE %.0f-%.0f keV",BATSE4,BATSE5);
  legL->AddEntry(Lct_BATSE4,legname);
  sprintf(legname, "GBM %.0f-%.0f keV",GBM1,GBM2);
  legL->AddEntry(Lct_GBM,legname);

  sprintf(legname, "SWIFT (BAT) %.0f-%.0f keV",SWIFT1,SWIFT2);
  legL1->AddEntry(Lct_SWIFT,legname);

  TLegend *legH = new TLegend(0.7,0.8,0.99,0.99);
  clc->cd(3);
  if(ExtractPhotons)
    {
      Lc->SetMarkerStyle(20);
      Lc->SetMarkerColor(2);
      Lc->SetLineColor(2);
      Lc->SetFillColor(2);
      Lc->SetStats(1);
      Lc->Draw();
      Lct_EXT->Draw("lsame");
      Lct_LAT->Draw("lsame");
      if(enph!=LAT1)    legH->AddEntry(Lct_EXT,"Extracted Photons (Theory)");
      legH->AddEntry(Lc,"Extracted Photons (Sampled)");
      std::cout<<" Photons Extracted : = "<<Lc->GetEntries()<<std::endl;
      std::cout<<" flux (exp)        : = "<<Counts->Integral(0,EBIN,"width")<<" keV"<<std::endl;
      std::cout<<" flux (exp)        : = "<<Counts->Integral(0,EBIN,"width")*1.0e-3/erg2meV/GenerationArea<<" erg/cm^2"<<std::endl;
    } 
  else
    {
      Lct_LAT->Draw("l");
    }
  
  legH->AddEntry(Lct_LAT,"LAT (tot)");
  clc->cd(1);
  legL1->Draw();
  clc->cd(2);
  legL->Draw();
  clc->cd(3);
  legH->Draw();
  
  std::cout<<" Nph LAT ("<<LAT1<<","<<LAT2<<")  = "<<sp->Integral_T(Lct_LAT,TMIN,TMAX)<<std::endl;
  if(ExtractPhotons) 
    {
      std::cout<<" Nph EXT ("<<enph<<","<<EMAX<<")  = "<<sp->Integral_T(Lct_EXT,TMIN,TMAX)<<std::endl;
    }

  std::cout<<" -------------------------------------------------- "<<std::endl;
  

  double fBATSE1  = sp->GetFluence(BATSE1,BATSE2);
  double fBATSE2  = sp->GetFluence(BATSE2,BATSE3);
  double fBATSE3  = sp->GetFluence(BATSE3,BATSE4);
  double fBATSE4  = sp->GetFluence(BATSE4,BATSE5);
  double fBATSE23 = sp->GetFluence(BATSE2,BATSE4);


  double fBATSET = sp->GetFluence(BATSE1,BATSE5);
  double fBAT    = sp->GetFluence(SWIFT1,SWIFT2);
  double fLAT   = sp->GetFluence(LAT1,LAT2);
  double fGBM   = sp->GetFluence(GBM1,GBM2);
  double fEXP   = sp->GetFluence(enph,EMAX);
  double fTOT   = sp->GetFluence();
  double fPeakFlux = sp->GetPeakFlux(BATSE2,BATSE4);

  double Ep = e2Ne->GetBinCenter(e2Ne->GetMaximumBin());
  std::cout<<"* Theoretical values:  *****************************"<<std::endl;
  std::cout<<" T90 = "<<T90<<" T100 = "<<T100<<" Epeak = "<<Ep<<std::endl;

  std::cout<<" GRB   fluence ("<< EMIN <<","<< EMAX <<") = "<<fTOT<<" erg/cm^2, "<<fTOT/T100<<" erg/cm^2/s "<<std::endl;
  std::cout<<" BATSE fluence (ch1) ("<<BATSE1<<","<<BATSE2<<") = "<<fBATSE1<<" erg/cm^2, "<<fBATSE1/T100<<" erg/cm^2/s "<<std::endl;
  std::cout<<" BATSE fluence (ch2) ("<<BATSE2<<","<<BATSE3<<") = "<<fBATSE2<<" erg/cm^2, "<<fBATSE2/T100<<" erg/cm^2/s "<<std::endl;
  std::cout<<" BATSE fluence (ch3) ("<<BATSE3<<","<<BATSE4<<") = "<<fBATSE3<<" erg/cm^2, "<<fBATSE3/T100<<" erg/cm^2/s "<<std::endl;
  std::cout<<" BATSE fluence (ch4) ("<<BATSE4<<","<<BATSE5<<") = "<<fBATSE4<<" erg/cm^2, "<<fBATSE4/T100<<" erg/cm^2/s "<<std::endl;
  std::cout<<" BATSE fluence (tot) ("<<BATSE1<<","<<BATSE5<<") = "<<fBATSET<<" erg/cm^2, "<<fBATSET/T100<<" erg/cm^2/s "<<std::endl;
  std::cout<<" ------------------------------ normalization: ------------------------------"<<std::endl;
  std::cout<<" BATSE Peak Flux ("<<BATSE2<<","<<BATSE4<<") = "<<fPeakFlux<<" ph/cm^2/s"<<std::endl;
  std::cout<<" BATSE Fluence   ("<<BATSE2<<","<<BATSE4<<") = "<<fBATSE23<<" erg/cm^2, "<<(fBATSE23)/T100<<" erg/cm^2/s "<<std::endl;
  std::cout<<" ------------------------------ normalization: ------------------------------"<<std::endl;
  std::cout<<" BAT   fluence ("<< SWIFT1 <<","<< SWIFT2 <<") = "<<fBAT<<" erg/cm^2, erg/cm^2/s "<<fBAT/T100<<std::endl;
  std::cout<<" GBM   fluence ("<< GBM1 <<","<< GBM2 <<") = "<<fGBM<<" erg/cm^2, erg/cm^2/s "<<fGBM/T100<<std::endl;
  std::cout<<" LAT   fluence ("<< LAT1 <<","<< LAT2 <<") = "<<fLAT<<" erg/cm^2, erg/cm^2/s "<<fLAT/T100<<std::endl;
  if(ExtractPhotons) std::cout<<" EXP   flux ("<< enph <<","<< EMAX <<") = "<<fEXP<<" erg/cm^2"<<std::endl; 
  animatedCanvas->cd();
  name2+=extension;
  if(savePlots)
    {
      
      cNv->Print(name2+"_Nv.eps");
      csp->Print(name2+"_sp.eps");
      clc->Print(name2+"_lc.eps");
    }

}

void MakeGRB(int NGRB=1, double enph=0, bool gbm=false)
{
  std::cout<<" ****** GRB and ROOT test ****** "<<std::endl;
  
  std::string path = ::getenv("GRBOBSROOT");
  std::string paramFile = path+"/src/test/GRBParam.txt";
  
  GRBobsParameters *params = new GRBobsParameters();  
  //////////////////////////////////////////////////
  params->ReadParametersFromFile(paramFile,NGRB);
  params->PrintParameters();
  GRBobsSim* m_grb = new GRBobsSim(params);

  double z = params->GetRedshift();
  double Eco = params->GetCutOffEnergy();
  if(ExtraComponent) 
    {
      TH2D *h = m_grb->MakeGRB_ExtraComponent(10000,nlat);
      m_grb->CutOff(h,Eco);
      m_grb->SaveNvEC();
    }
  else
    {
      TH2D *h = m_grb->MakeGRB();
      m_grb->CutOff(h,Eco);
      m_grb->SaveNv();
    }
  char GRBname[100];
  sprintf(GRBname,"%d",(int) params->GetGRBNumber());
  if (gbm)  m_grb->GetGBMFlux(GRBname);

  delete m_grb;
  char name[100];
  sprintf(name,"grbobs_%d.root",(int) params->GetGRBNumber());
  if(ExtraComponent)  sprintf(name,"grbobs_%d_EC.root",(int) params->GetGRBNumber());

  extension=params->GetGRBNumber();

  delete params; //??


  PlotGRB(enph,z,name);
  
 }

//////////////////////////////////////////////////
void ScanParameters(int Ngrb)
{
  
  std::string path = ::getenv("GRBOBSROOT");
  std::string paramFile = path+"/src/test/GRBParam.txt";  
  GRBobsParameters *params = new GRBobsParameters();  

  //////////////////////////////////////////////////
  // OUTPUT
  double T90;
  double fPeakFlux;
  double fBATSE1,fBATSE2,fBATSE3,fBATSE4;
  double fBATSE, fLAT,fGBM,fEXP,fTOT;
  double nLAT,nLAT30,nLAT100,nLAT1000;
  //  double nBATSE,nGBM,nEXP;
  double alpha,beta,Ep;
  //////////////////////////////////////////////////
  TTree *GRBTree = new TTree("GRBTree","GRBOBS Catalogue");
  GRBTree->Branch("T90",&T90,"T90/D");
  
  GRBTree->Branch("LogfBATSE",&fBATSE,"fBATSE/D");  
  GRBTree->Branch("Logf1BATSE",&fBATSE1,"fBATSE1/D");
  GRBTree->Branch("Logf2BATSE",&fBATSE2,"fBATSE2/D");
  GRBTree->Branch("Logf3BATSE",&fBATSE3,"fBATSE3/D");
  GRBTree->Branch("Logf4BATSE",&fBATSE4,"fBATSE4/D");

  GRBTree->Branch("LogPeakFlux",&fPeakFlux,"fPeakFlux/D");

  GRBTree->Branch("LogfLAT",&fLAT,"LogfLAT/D");
  GRBTree->Branch("LogfGBM",&fGBM,"LogfGBM/D");
  GRBTree->Branch("LogfEXP",&fEXP,"LogfEXP/D");
  GRBTree->Branch("LogfTOT",&fTOT,"LogfTOT/D");
  
  //  GRBTree->Branch("LognBATSE",&nBATSE,"LognBATSE/D");
  //  GRBTree->Branch("LognGBM",&nGBM,"LognGBM/D");

  GRBTree->Branch("LognLAT",&nLAT,"LognLAT/D");
  GRBTree->Branch("LognLAT30",&nLAT30,"LognLAT30/D");
  GRBTree->Branch("LognLAT100",&nLAT100,"LognLAT100/D");
  GRBTree->Branch("LognLAT1000",&nLAT1000,"LognLAT1000/D");

  //////////////////////////////////////////////////
  GRBTree->Branch("LogEp",&Ep,"LogEp/D");
  GRBTree->Branch("alpha",&alpha,"alpha/D");
  GRBTree->Branch("beta",&beta,"beta/D");
  
  for(long grbn=1;grbn<=Ngrb;grbn++)
    {
      std::cout<<" GRBOBS number : "<<grbn<<" / "<<Ngrb<<std::endl;
      params->ReadParametersFromFile(paramFile,grbn);
      GRBobsSim m_grb(params);
      TH2D *Nv = m_grb.MakeGRB();
      
      TMIN = Nv->GetXaxis()->GetXmin();
      TMAX = Nv->GetXaxis()->GetXmax();
      EMIN = Nv->GetYaxis()->GetXmin();
      EMAX = Nv->GetYaxis()->GetXmax();
      EBIN = Nv->GetYaxis()->GetNbins();
      double dt = Nv->GetXaxis()->GetBinWidth(1);
      
      // e2Ne = e * Ne * de: [keV/s/m^2]
      gDirectory->Delete("e2Ne");
      gDirectory->Delete("LAT");
      Nv->ProjectionY("e2Ne");
      Nv->ProjectionY("LAT");

      TH1D *e2Ne = (TH1D*) gDirectory->Get("e2Ne");
      TH1D *LAT = (TH1D*) gDirectory->Get("LAT");
      Ep=0.0;
      for(int i=0; i < EBIN; i++)
	{
	  double ne = e2Ne->GetBinContent(i+1);
	  double de = e2Ne->GetBinWidth(i+1);//GetBinCenter(i+1);
	  double en = e2Ne->GetBinCenter(i+1);

	  double Aeff_on = AeffLAT(en*1e-3,1.0); // cm^2
	  LAT->SetBinContent(i+1, dt * de * ne * Aeff_on); // ph
	  e2Ne->SetBinContent(i+1,en * en * ne);
	}
      Ep = log10(e2Ne->GetBinCenter(e2Ne->GetMaximumBin()));
      alpha =  params->GetAlpha();
      beta =   params->GetBeta();
      //////////////////////////////////////////////////
      SpectObj *sp = new SpectObj(Nv,0, params->GetRedshift());              //ph
      sp->SetAreaDetector(GenerationArea); //like observation sim
      //////////////////////////////////////////////////
  
      TH1D *Lct_TOT   = sp->Integral_E(EMIN,EMAX);  // ph (over generation area)
      TH1D *Lct_BATSE = sp->Integral_E(BATSE1,BATSE5);  // ph
      TH1D *Lct_GBM   = sp->Integral_E(GBM1,GBM2);  // ph

      TH1D *Lct_LAT30   = sp->Integral_E(30000.,LAT2);  // ph
      TH1D *Lct_LAT100   = sp->Integral_E(100000.,LAT2);  // ph
      TH1D *Lct_LAT1000   = sp->Integral_E(1000000.,LAT2);  // ph
      //      TH1D *Lct_EXT   = sp->Integral_E(enph,EMAX);  // ph
      
      T90    = log10(sp->GetT90());
      fTOT   = log10(sp->GetFluence());
      fBATSE = log10(sp->GetFluence(BATSE1,BATSE5));
      fBATSE1 = log10(sp->GetFluence(BATSE1,BATSE2));
      fBATSE2 = log10(sp->GetFluence(BATSE2,BATSE3));
      fBATSE3 = log10(sp->GetFluence(BATSE3,BATSE4));
      fBATSE4 = log10( sp->GetFluence(BATSE4,BATSE5));

      fPeakFlux = log10(sp->GetPeakFlux(BATSE2,BATSE4));

      fGBM   = log10(sp->GetFluence(GBM1,GBM2));
      fLAT   = log10(sp->GetFluence(LAT1,LAT2));
      fEXP   = log10(sp->GetFluence(enph,emax));
      
      nLAT   = log10(LAT->Integral());
      //      nBATSE = log10(sp->Integral_T(Lct_BATSE,0.0,TMAX));
      //      nGBM   = log10(sp->Integral_T(Lct_GBM,0.0,TMAX));

      nLAT30   = log10(sp->Integral_T(Lct_LAT30,TMIN,TMAX));
      nLAT100   = log10(sp->Integral_T(Lct_LAT100,TMIN,TMAX));
      nLAT1000   = log10(sp->Integral_T(Lct_LAT1000,TMIN,TMAX));

      double Ep = e2Ne->GetBinCenter(e2Ne->GetMaximumBin());
      std::cout<<"* Theoretical values:  *****************************"<<std::endl;
#ifdef WIN32
      std::cout<<" (not available in windows)" << std::endl;
#else
      std::cout<<" T90 = "<<pow(10.,T90)<<" Epeak = "<<Ep<<std::endl;
      std::cout<<" log TOT   flux ("<< EMIN <<","<< EMAX <<") = "<<fTOT<<" erg/cm^2"<<std::endl;
      std::cout<<" BATSE flux (ch1) ("<<BATSE1<<","<<BATSE2<<") = "<<pow(10.,fBATSE1)<<" erg/cm^2"<<std::endl;
      std::cout<<" BATSE flux (ch2) ("<<BATSE2<<","<<BATSE3<<") = "<<pow(10.,fBATSE2)<<" erg/cm^2"<<std::endl;
      std::cout<<" BATSE flux (ch3) ("<<BATSE3<<","<<BATSE4<<") = "<<pow(10.,fBATSE3)<<" erg/cm^2"<<std::endl;
      std::cout<<" BATSE flux (ch4) ("<<BATSE4<<","<<BATSE5<<") = "<<pow(10.,fBATSE4)<<" erg/cm^2"<<std::endl;
      std::cout<<" BATSE flux (tot) ("<<BATSE1<<","<<BATSE5<<") = "<<pow(10.,fBATSE)<<" erg/cm^2"<<std::endl;
      std::cout<<" BATSE Peakflux ("<<BATSE2<<","<<BATSE4<<") = "<<pow(10.,fPeakFlux)<<" ph/cm^2/s"<<std::endl;
      std::cout<<" GBM   flux ("<< GBM1 <<","<< GBM2 <<") = "<<pow(10.,fGBM)<<" erg/cm^2"<<std::endl;
      std::cout<<"  LAT  flux ("<< LAT1 <<","<< LAT2 <<") = "<<pow(10.,fLAT)<<" erg/cm^2"<<std::endl;
      //      std::cout<<"  Nph BATSE  ("<<BATSE1<<","<<BATSE5<<") = "<<pow(10.,nBATSE)<<std::endl;
      //      std::cout<<"  Nph GBM    ("<<GBM1<<","<<GBM2<<")  = "<<pow(10.,nGBM)<<std::endl;
      std::cout<<"  Nph LAT    ("<<EMIN<<","<<EMAX<<")  = "<<pow(10.,nLAT)<<std::endl;
      std::cout<<"  Nph LAT30    ("<<30000<<","<<LAT2<<")  = "<<pow(10.,nLAT30)<<std::endl;
      std::cout<<"  Nph LAT100    ("<<100000<<","<<LAT2<<")  = "<<pow(10.,nLAT100)<<std::endl;
      std::cout<<"  Nph LAT1000    ("<<1000000<<","<<LAT2<<")  = "<<pow(10.,nLAT1000)<<std::endl;
#endif
      //	  if(output) delete sp;
      GRBTree->Fill();
      if(grbn%10==0)
	{
	  TFile aFile("GRBobsCatalogueFile.root","RECREATE");
	  GRBTree->Write();
	}
      delete sp;
    }
  TFile aFile("GRBobsCatalogueFile.root","RECREATE");
  GRBTree->Write();
  delete params;
}


int main(int argc, char** argv)
{
  
  std::string arg_name("");
  int current_arg = 1;
  double enph=0.0;
  int ngrb=1;
  int ngrbs=0;
  bool gbm = false;
  
  while(current_arg < argc)
    {
      arg_name = argv[current_arg];
      
      if("-extract"==arg_name)
	{
	  enph  = atof(argv[++current_arg]);
	}
      else if("-grb"==arg_name)
	{
	  ngrb=atoi(argv[++current_arg]);
	}
      else if("-gbm"==arg_name)
	{
	  gbm=true;
	}
      else if ("-EC" ==arg_name)
	{
	  ExtraComponent=true;
	  nlat=atoi(argv[++current_arg]);
	}
      else if("-scan"==arg_name)
	{
	  ngrbs=atoi(argv[++current_arg]);
	}
      else if("-save"==arg_name)
	{
	  savePlots=true;
	  frame=atoi(argv[++current_arg]);
	}
      else if("-band"==arg_name)
	{
	  bandFit=true;
	}
      else if("-powerLaw"==arg_name)
	{
	  powerlawFit=true;
	}
      else if("-scaled"==arg_name)
	{
	  scaled=true;
	}
      else if("-noise"==arg_name)
	{
	  noise=true;
	  noiseL = TMath::Max(0.0,atof(argv[++current_arg]));
	}
      else if("-movie"==arg_name)
	{
	  movie=true;
	}
      
      current_arg++;
    }
  
  if(ngrbs>0)
    {
      ScanParameters(ngrbs);
    }
  else
    {
      TApplication theApp("App",0,0);
      MakeGRB(ngrb,enph,gbm);
      theApp.Run();
    }
}

void help()
{
  std::cout<<"--------------------------------------------------"<<std::endl;
  std::cout<<"   GRB ROOT test program author: Nicola Omodei    "<<std::endl;
  std::cout<<"   Option are: "<<std::endl;
  std::cout<<"  -extract [enph] : etxtract photons above enph"<<std::endl;
  std::cout<<"  -grb [N] processes the N grb in the $GRBROOT/src/test/GRBParam.txt file "<<std::endl;
  std::cout<<"  -band    Fit the GBM spectrum with as Band function, as afunction of the time"<<std::endl;
  std::cout<<"--------------------------------------------------"<<std::endl;

}

