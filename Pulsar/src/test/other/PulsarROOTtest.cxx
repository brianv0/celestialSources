#include <iostream>
#include <vector>

#include "TROOT.h"
#include "TApplication.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "SpectObj/SpectObj.h"
#include "Pulsar/PulsarSim.h"
#include "Pulsar/PulsarConstants.h"

using namespace cst;
double EMIN, EMAX, TMIN, TMAX, DT;
int    TBIN, EBIN;
//////////////////////////////////////////////////

TH2D *Load(char name[100]="pulsar.root")
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
  Nv->SetZTitle("N_{v} [ph/m^2/s/keV]");
  Nv->GetXaxis()->SetTitleOffset(1.5);
  Nv->GetYaxis()->SetTitleOffset(1.5);
  Nv->GetZaxis()->SetTitleOffset(1.2);
  Nv->GetXaxis()->CenterTitle();
  Nv->GetYaxis()->CenterTitle();
  Nv->GetZaxis()->CenterTitle();
  
  return Nv;
}



void MakePulsar()
{
}//[in-charge]


void PlotPulsar(double enph = 0,char name[100]="pulsar.root")
{
  gStyle->SetCanvasColor(10);
  
  TCanvas *cNv = new TCanvas("cNv","Nv");
  cNv->SetLogy();
  cNv->SetLogz();
  
  TH2D *Nv = Load(name);

  // Ne =  [ph/keV/s/m²]
  gDirectory->Delete("Ne");
  Nv->ProjectionY("Ne");
  TH1D *Ne = (TH1D*) gDirectory->Get("Ne");
  //  Ne->Scale(DT);
  Ne->SetYTitle("Ne [ph/MeV/cm^2/s]");
  Ne->SetXTitle(" Energy[keV] ");
  // Fv = Ne * de: [ph/s/m^2]
  TH1D *Fv = (TH1D*) Ne->Clone();
  Fv->SetTitle("Fv");
  Fv->SetName("Fv");
  Fv->SetYTitle("Fv [MeV/MeV/cm^2/s]");
  
  // vFv = e * Ne * de: [keV/s/m^2]
  TH1D *vFv = (TH1D*) Ne->Clone();
  vFv->SetTitle("vFv");
  vFv->SetName("vFv");

  vFv->SetYTitle(" vFv ");
  vFv->SetStats(0);

  Ne->SetLineColor(2);
  Fv->SetLineColor(3);
  vFv->SetLineColor(4);

  Ne->SetLineWidth(2);
  Fv->SetLineWidth(2);
  vFv->SetLineWidth(2);
  

  for(int i=0; i < EBIN; i++)
    {
      double ne = Ne->GetBinContent(i+1);
      double de = Ne->GetBinWidth(i+1);//GetBinCenter(i+1);
      double e = Ne->GetBinCenter(i+1);
      Fv->SetBinContent(i+1,de*ne);
      vFv->SetBinContent(i+1,e*de*ne);
    }
  
  

  //////////////////////////////////////////////////
  
  Nv->Draw("surf");
  
  bool ExtractPhotons = true;
  if(enph<=0) ExtractPhotons = false;
  else if (enph<EMIN)  enph = EMIN;
  
  //////////////////////////////////////////////////
  SpectObj *sp = new SpectObj(Nv,1);              //ph/m² 

  //////////////////////////////////////////////////
  TCanvas *clc = new TCanvas("clc","clc",600,800);
  TCanvas *csp = new TCanvas("csp","csp");
 
  csp->SetLogx();
  csp->SetLogy();
  
  clc->Divide(1,2);
  TH1D *Lct_Pulsar = sp->Integral_E(EMIN,EMAX);  // ph/m^2
  TH1D *Lct_EGRET = sp->Integral_E(EGRET2,EGRET3);  // ph/m^2
 
  TH1D *Lct_LAT = sp->Integral_E(LAT1,LAT2);  // ph/m^2
  TH1D *Lct_EXT = sp->Integral_E(enph,EMAX);  // ph/m^2
  
  Lct_Pulsar->SetStats(0);
  Lct_LAT->SetStats(0);
  
  Lct_EGRET->SetTitle("flux [ph]");
  Lct_EGRET->SetXTitle("Time (s)");
  Lct_EGRET->SetYTitle("photons/m^2");
    
  Lct_Pulsar->SetName("flux Pulsar [ph]");
  Lct_EGRET->SetName("flux EGRET [ph]");
  Lct_LAT->SetName("flux LAT [ph]"); 
  Lct_EXT->SetName("flux EXT [ph]"); 
 
  Lct_Pulsar->SetLineWidth(2);
  Lct_EGRET->SetLineWidth(2);
  Lct_LAT->SetLineWidth(2);
  Lct_EXT->SetLineWidth(2);
  
  Lct_Pulsar->SetLineColor(2);
  Lct_EGRET->SetLineColor(3);
  Lct_LAT->SetLineColor(4);
  Lct_EXT->SetLineColor(6);
 
  csp->cd();
  vFv->SetMinimum(.1);

  /* these lines are to be erased if we want ph-kev and not ph-kev-sec

  */
  double Period = 0.089;  
  Fv->Scale(.1);///Period);
  Ne->Scale(1e-4);///Period);
  //  vFv->Draw("al");
  Fv->Draw("al");
  Ne->Draw("samel");

  std::cout << " --- \nNe @ 100 MeV " << Ne->GetBinContent(Ne->FindBin(EGRET2)) << " ph/cm2/Mev/s \n---" << std::endl;

  //  gDirectory->Delete("band");
  
  TLegend *leg = new TLegend(0.11,0.12,0.37,0.25);
  leg->SetFillStyle(0);
  leg->AddEntry(Ne," Ne  [ph/MeV/cm^2/s] ","lp");
  leg->AddEntry(Fv," Fv  [ph/cm^2/s] ","lp");
  leg->AddEntry(vFv," vFv [MeV/cm^2/s] ","lp");
  leg->Draw();
  
  TH1D *Counts   = sp->GetSpectrum(); //ph
  Counts->SetName("Counts");
  TH1D *Lc = new TH1D("Lc","Lc",100,0,10000*TMAX);//sp->GetTimes();  
  Lc->SetName("Counts[ph]");
  
  if(ExtractPhotons)
    {  
      double time =0.0;
      double Interval;
      double energy;
      double flux;
      int i = 0;
      time = sp->interval(0.0,enph); // s
      while(time < 10000*TMAX && time>=0)
	{
	  Interval = sp->interval(time,enph); // s  
	  energy   = sp->energy(time,enph);   // keV
	  flux     = sp->flux(time,enph);     // ph/s
	  
	  Counts->Fill(energy); // ph
	  Lc->Fill(time);       // ph
	  
	  std::cout<<" Time (s)  = "<<time
		   <<" Flux (ph/s/m^2)  = "<<flux
		   <<" energy (keV) = "<<energy
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
      Counts->Draw("E1same");
      Counts->SetStats(1);
    }
  //////////////////////////////////////////////////      
  //clc->cd();
  double LowEnergyMax  = 1.05*Lct_Pulsar->GetMaximum();
  double HighEnergyMax = Lct_LAT->GetMaximum();
  double tmax = Lct_Pulsar->GetXaxis()->GetXmax();
  double tmin = Lct_Pulsar->GetXaxis()->GetXmin();
  
  clc->cd(1);
  Lct_Pulsar->Draw("l");
  Lct_EGRET->Draw("lsame");

  clc->cd(2);
  
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
      std::cout<<" Photons Extracted : = "<<Lc->GetEntries()<<std::endl;
      std::cout<<" flux (exp)        : = "<<Counts->Integral(0,EBIN,"width")<<" keV/m^2"<<std::endl;
      std::cout<<" flux (exp)        : = "<<Counts->Integral(0,EBIN,"width")*1.0e-7/erg2meV<<" erg/cm^2"<<std::endl;
      
    } 
  else
    {
      Lct_LAT->Draw("l");
    }



  std::cout<<" Nph Pulsar ("<<EMIN<<","<<EMAX<<")  = "<<sp->Integral_T(Lct_Pulsar,0.0,TMAX)<<std::endl;
  std::cout<<" Nph EGRET ("<<EGRET2<<","<<EGRET3<<")  = "<<sp->Integral_T(Lct_EGRET,0.0,TMAX)<<std::endl;
  std::cout<<" Nph LAT ("<<LAT1<<","<<LAT2<<")  = "<<sp->Integral_T(Lct_LAT,0.0,TMAX)<<std::endl;
  std::cout<<" Nph EGRET ("<<EGRET2<<","<<EGRET3<<")  = "<<sp->Integral_T(Lct_EGRET,0.0,TMAX)*1e-4/TMAX<< " ph/cm2/s " << std::endl;


  if(ExtractPhotons) std::cout<<" Nph EXT ("<<enph<<","<<EMAX<<")  = "<<sp->Integral_T(Lct_EXT,0.0,TMAX)<<std::endl;
  std::cout<<" -------------------------------------------------- "<<std::endl;
  
  int iEGRET2 = Ne->FindBin(EGRET2);
  int iEGRET3 = Ne->FindBin(EGRET3);
  int iLAT1 = Ne->FindBin(LAT1);
  int iLAT2 = Ne->FindBin(LAT2);
  int iEXP = Ne->FindBin(enph);
  
  std::cout<<" Flux[ erg/cm^2] Pulsar ("<<EMIN<<","<<EMAX<<")  = "<<Fv->Integral(0,EBIN,"width")*1.0e-7/erg2meV<<" erg/cm^2"<<std::endl;
  std::cout<<" Flux[ erg/cm^2] EGRET ("<<EGRET2<<","<<EGRET3<<")  = "<<Fv->Integral(iEGRET2,iEGRET3,"width")*1.0e-7/erg2meV<<" erg/cm^2"<<std::endl;
  std::cout<<" Flux[ erg/cm^2] LAT ("<<LAT1<<","<<LAT2<<")  = "<<Fv->Integral(iLAT1,iLAT2,"width")*1.0e-7/erg2meV<<" erg/cm^2"<<std::endl;  
  if(ExtractPhotons) std::cout<<" Flux[ erg/cm^2] EXT ("<<enph<<","<<EMAX<<")  = "<<Fv->Integral(iEXP,EBIN,"width")*1.0e-7/erg2meV<<" erg/cm^2"<<std::endl;
  
}


int main(int argc, char** argv)
{
  
  std::cout<<" ****** Pulsar and ROOT test ****** "<<std::endl;
  std::string arg_name("");
  int current_arg = 1;
  double enph=0.0;
  while(current_arg < argc)
    {
      arg_name = argv[current_arg];
      
      if("-extract"==arg_name)
	{
	  enph  = atof(argv[++current_arg]);
	}
      current_arg++;
    }
  TApplication theApp("App",0,0);

  /*
  double Period  = 0.089; // s
  double Fluence = 9e-6; // ph/cm2/s
  int npeaks = 2;
  double ppar1 = 1e6;
  double ppar2 = 8e6;
  double ppar3 = -1.62;
  double ppar4 = 1.7;
  */

  
  //Crab Polar Cap vs Outer Gap
  double Period  = 0.033; // s
  double Fluence = 2.3e-6; // ph/cm2/s
  int npeaks = 2;
  double ppar1 = 1e6;
  double ppar2 = 30e6;
  double ppar3 = -1.9;
  double ppar4 = 0.29; // 1.0 for outer
  

  /*
  //Vela Polar Cap vs Outer 
  double Period  = 0.089; // s
  double Fluence = 9e-6; // ph/cm2/s
  int npeaks = 2;
  double ppar1 = 1e6;
  double ppar2 = 8e6;
  double ppar3 = -1.62;
  double ppar4 = 1.7; // for outer
  //double ppar4 = 1.7; // for polar
  */ 
  /*
  //Geminga Polar Cap vs Outer 
  double Period  = 0.237; // s
  double Fluence = 366e-6; // ph/cm2/s
  int npeaks = 2;
  double ppar1 = 1e6;
  double ppar2 = 4.8e6;
  double ppar3 = -1.56;
  double ppar4 = 1.0; // 1.0 for outer
 */



  std::cout << " --> Pulsar Sim init " << Fluence << " " << Period << std::endl; 
  PulsarSim* m_pulsar = new PulsarSim(Fluence,Period, npeaks);
  m_pulsar->SaveNv((TH2D*)m_pulsar->PSRPhenom(ppar1,ppar2,ppar3,ppar4));
  char name[100];
  sprintf(name,"pulsar.root");
  PlotPulsar(enph,name);
  
  theApp.Run();
}
