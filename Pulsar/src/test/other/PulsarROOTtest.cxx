
#include <iostream>
#include <vector>
#include <fstream>
#include <time.h>

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
#include <flux/EventSource.h>

using namespace cst;
double EMIN, EMAX, TMIN, TMAX, DT;
int    TBIN, EBIN;
double AreaDetector = 0.0; //m2
const double nLoops = 970787;//=1day of Vela;

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
  DT=TMAX/TBIN;
  std::cout<<"****  Nv Matrix Loaded ! "<<std::endl; 
  std::cout<<"      Bin t = "<<TBIN<<": TMIN = "<<TMIN<<" TMAX = "<<TMAX<<std::endl; 
  std::cout<<"      Bin e = "<<EBIN<<": EMIN = "<<EMIN<<" EMAX = "<<EMAX<<std::endl; 
  
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

  std::cout << "\n**** Plotting Pulsar Simulator Results..." << std::endl;
  gStyle->SetCanvasColor(10);
  
  TCanvas *cNv = new TCanvas("cNv","Nv");
  cNv->SetLogy();
  cNv->SetLogz();
  
  TH2D *Nv = Load(name); //Nv = ph/m2/s/keV

  // Ne =  [ph/m²/s/keV]
  gDirectory->Delete("Ne");
  Nv->ProjectionY("Ne");


  TH1D *Ne = (TH1D*) gDirectory->Get("Ne");

  Ne->SetYTitle("Ne [ph/m^2/s/keV]");
  Ne->SetXTitle(" Energy[keV] ");
  // Fv = Ne * de: [ph/m^2/s]
  TH1D *Fv = (TH1D*) Ne->Clone();
  Fv->SetTitle("Fv");
  Fv->SetName("Fv");
  Fv->SetYTitle("Fv [ph/m^2/s]");
  
  // vFv = e * Ne * de: [keV/m^2/s]
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
  
  Nv->SetMinimum(1e-20);
  Nv->Draw("surf");
  
  bool ExtractPhotons = true;
  if(enph<=0) ExtractPhotons = false;
  else if (enph<EMIN)  enph = EMIN;
  
  //////////////////////////////////////////////////
  SpectObj *sp = new SpectObj(Nv,1);              //ph 
  AreaDetector = 1.21;//EventSource::totalArea();
  sp->SetAreaDetector(1.21);//EventSource::totalArea());
  std::cout << "**  Effective Area set to : " << sp->GetAreaDetector() << " m2 " << std::endl; 
  //////////////////////////////////////////////////
  TCanvas *clc = new TCanvas("LightCurves","LightCurves",600,800);
  clc->Divide(1,2);

  TCanvas *csp = new TCanvas("csp","csp");
  csp->SetLogx();
  csp->SetLogy();

  //Init lightCurves

  //Total band
  TH1D *Lct_Pulsar = sp->Integral_E(EMIN,EMAX);  // ph
  Lct_Pulsar->Scale((1e-4)/(DT*AreaDetector)); //ph/cm2/s
  //Normalisation band
  TH1D *Lct_NORM = sp->Integral_E(EnNormMin,EnNormMax);  // ph
  Lct_NORM->Scale((1e-4)/(DT*AreaDetector)); //ph/cm2/s
  //EGRET band
  TH1D *Lct_EGRET = sp->Integral_E(EGRET1,EGRET3);  // ph
  Lct_EGRET->Scale((1e-4)/(DT*AreaDetector)); //ph/cm2/s
  //LAT band
  TH1D *Lct_LAT = sp->Integral_E(LAT1,LAT2);  // ph
  Lct_LAT->Scale((1e-4)/(DT*AreaDetector)); //ph/cm2/s
  //Extracted photons band
  TH1D *Lct_EXT = sp->Integral_E(enph,EMAX);  // ph

  //  std::cout << "### " << Lct_NORM->GetBinContent(12) << std::endl;

  Lct_Pulsar->SetStats(0);
  Lct_LAT->SetStats(0);
  Lct_EGRET->SetStats(0);
  Lct_NORM->SetStats(0);
  

  Lct_Pulsar->SetTitle("Flux [ph/cm2/s]");
  Lct_Pulsar->SetXTitle("Time (s)");
  Lct_Pulsar->SetYTitle("ph/cm^2/s");
    
  Lct_Pulsar->SetName("Total Flux  [ph/cm2/s]");
  Lct_EGRET->SetName("Flux EGRET band [ph]");
  Lct_LAT->SetName("Flux LAT band [ph]"); 
 
  Lct_Pulsar->SetLineWidth(2); 
  Lct_EGRET->SetLineWidth(2);
  Lct_LAT->SetLineWidth(2);
  
  Lct_Pulsar->SetLineColor(2);
  Lct_EGRET->SetLineColor(7);
  Lct_LAT->SetLineColor(4);
 

  csp->cd();
  
  vFv->Scale(1e-6); //Scale in order to plot the 3 histos on the same canvas

  Ne->SetMinimum(1e-10);
  Fv->SetMinimum(1e-10);
  vFv->SetMinimum(1e-10);

  Fv->Draw("al");
  vFv->Draw("alsame");
  Ne->Draw("samel");
  
  std::cout << "****  Ne @ 100 MeV " << Ne->GetBinContent(Ne->FindBin(EGRET2)) << " ph/cm2/kev/s " << std::endl;

  //  gDirectory->Delete("band");
  
  TLegend *leg = new TLegend(0.11,0.12,0.37,0.25);
  leg->SetFillStyle(0);
  leg->AddEntry(Ne," Ne  [ph/keV/m^2/s] ","lp");
  leg->AddEntry(Fv," Fv  [ph/m^2/s] ","lp");
  leg->AddEntry(vFv," vFv [keV/m^2/s] - Scaled by 1e6 -","lp");
  leg->Draw();
  
  TH1D *Counts   = sp->GetSpectrum(); //ph
  Counts->SetName("Counts");
  TH1D *Lc = new TH1D("Lc","Lc",TBIN,0,TMAX);//sp->GetTimes();  
  Lc->SetName("Counts[ph]");
  Lc->SetMarkerStyle(20);
  Lc->SetMarkerSize(0.5);
  Lc->SetMarkerColor(1);
  Lc->SetLineColor(1);
  Lc->SetStats(0);  


  if(ExtractPhotons)
    {  
      double time =0.0;
      double Interval;
      double energy;
      double flux;
      int i = 0;
      time = sp->interval(0.0,enph); // s
      std::cout << "Extracting photon above " << enph << std::endl;



      while(time < nLoops*TMAX && time>=0)
	{
	  Interval = sp->interval(time,enph); // s  
	  energy   = sp->energy(time,enph);   // keV
	  flux     = sp->flux(time,enph);     // ph/s
	  
	  Counts->Fill(energy); // ph
	  Lc->Fill(time - (TMAX*int(time/TMAX)));       // ph
	  time+=Interval;
	  /*
	    std::cout << " Phase " << (time - (TMAX*int(time/TMAX)))/TMAX << std::endl;
	    std::cout<<" Time (s)  = "<<time
	       << " within period = " << time - (TMAX*int(time/TMAX))
	       <<" Flux (ph/s/cm^2)  = "<<flux*1e-4
	       <<" energy (keV) = "<<energy
	       <<" Interval (s) = "<<Interval<<std::endl;
	  */ 
	
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
      //Counts->SetMaximum(1.5*Ne->GetMaximum());

      Ne->Scale(nLoops*DT*AreaDetector);
      Fv->Scale(nLoops*DT*AreaDetector);//*DT*6e4);
      vFv->Scale(nLoops*DT*AreaDetector);
      
      Ne->SetMinimum(1e-10*DT*AreaDetector);
      Fv->SetMinimum(1e-10*DT*AreaDetector);
      vFv->SetMinimum(1e-10*DT*AreaDetector);

      Fv->Draw("al");
      vFv->Draw("alsame");
      Ne->Draw("samel");

      Counts->SetStats(1);
      Counts->Draw("E1same");
      Counts->SetStats(1);
    }

  TLegend *lcleg = new TLegend(0.30,0.89,0.88,1.0);
  lcleg->SetFillStyle(0);
  lcleg->SetTextSize(0.03);
  lcleg->SetFillColor(10);
  lcleg->AddEntry(Lct_Pulsar,"Total flux in ph/cm2/s (30MeV-400GeV)","lp");
  lcleg->AddEntry(Lct_EGRET,"EGRET flux in ph/cm2/s (30Mev-30GeV)","lp");
  lcleg->AddEntry(Lct_LAT,"LAT flux in ph/cm2/s (30 MeV-300GeV)","lp");
  clc->cd(1);
  Lct_Pulsar->Draw();
  Lct_EGRET->Draw("same");
  
  clc->cd(2);

  std::cout << "****-------------------------------------------------" <<std::endl;
  std::cout << "Fluxes summary: (Eff.Area=" << AreaDetector << " m2 )" << std::endl;
  std::cout << "       Normalisation:   (" << EnNormMin <<","<< EnNormMax <<")  = "
	    << Lct_NORM->Integral(0,TBIN)*(DT/TMAX)  << " ph/cm2/s" << std::endl;
  std::cout << "       Total (" << EMIN<<","<<EMAX<<")  = "  
	    << Lct_Pulsar->Integral(0,TBIN)*(DT/TMAX) << " ph/cm2/s" << std::endl;
  std::cout << "       EGRET band  (" << EGRET1 <<","<<EGRET3<<")  = "
	    << Lct_EGRET->Integral(0,TBIN)*(DT/TMAX) << " ph/cm2/s" << std::endl;

  if(ExtractPhotons)
    {
      std::cout << "       LAT band (" << LAT1 <<","<<LAT2<<")  = " 
		<< Lct_LAT->Integral(0,TBIN)*(DT/TMAX) << " ph/cm2/s" << std::endl;


  
      int en2 = Fv->GetXaxis()->FindBin(EnNormMin);
      int en3 = Fv->GetXaxis()->FindBin(EnNormMax);
      std::cout<<" Photons Extracted : = "<<Lc->GetEntries()<<std::endl;
      //std::cout<<" flux (exp)        : = "<<Counts->Integral(en2,en3)/(1e4*nLoops*AreaDetector)<<" ph/m^2"<<std::endl;
      //std::cout<<" flux (exp)        : = "<<Counts->Integral(0,EBIN,"width")*1.0e-7/erg2meV<<" erg/cm^2"<<std::endl;
      Lct_LAT->Scale(nLoops*DT*6e4);

      lcleg->AddEntry(Lc,"LAT extracted photons in ph/cm2/s (30 MeV-300GeV)","lp");
      
      Lct_EXT->SetLineColor(8);
      Lct_EXT->Scale(nLoops);
      Lct_EXT->SetLineWidth(2);

      Lct_EXT->Draw(); 
      Lc->SetStats(0);
      Lct_EXT->SetStats(0);

      Lc->Draw("epsame");            

      std::cout << " Photon expecteed : " << Lct_EXT->Integral(0,TBIN) << " ph. " << std::endl;  

    } 
  else
    {
      std::cout << "       LAT band (" << LAT1 <<","<<LAT2<<")  = " 
		<< Lct_LAT->Integral(0,TBIN) << " ph/cm2/s" << std::endl;
      Lct_Pulsar->Draw("l");
      Lct_LAT->Draw("lsame");
    }

  lcleg->Draw("lsame");


  std::cout << "****-------------------------------------------------" <<std::endl;

 
  delete sp; 
}


int main(int argc, char** argv)
{
  
  std::cout<<"\n******** Pulsar and ROOT test ********"<<std::endl;
  std::string arg_name("");
  int current_arg = 1;
  double enph= 0.0;
  double enphmin = 3e4;
  double enphmax = 3e8;

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

  

  
  double Period  = 0.089; // s
  double flux = 1e-5; // ph/cm2/s
  int npeaks = 2;
  double ppar1 = 1e6;
  double ppar2 = 8e6;
  double ppar3 = -1.62;
  double ppar4 = 1.7;//1.7;
  int seed = 61443;


  PulsarSim* m_pulsar = new PulsarSim("PSRVELA",seed,flux,enphmin, enphmax, Period);
  m_pulsar->SaveNv((TH2D*)m_pulsar->PSRPhenom(double(npeaks), ppar1,ppar2,ppar3,ppar4));

  

  char name[100];
  sprintf(name,"PSRVELAroot.root");
  
  std::cout << "**  Photons simulated on a period of " << Period*nLoops << " s. " << std::endl;
  PlotPulsar(enph,name);  

  theApp.Run();
}
