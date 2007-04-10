#include <iostream>
#include <vector>

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

#include "facilities/Util.h"

//#include "../src/GRB/SpectObj.h"
#include "../../GRB/SpectObj.h"
#include "../../GRB/GRBSim.h"
#include "../../GRB/GRBConstants.h"

using namespace cst;
double EMIN, EMAX, TMIN, TMAX, DT;
int    TBIN, EBIN;
//////////////////////////////////////////////////

Parameters *params;

double Band(double *var, double *par)
{

  double a  = par[0];
  double b  = par[1];
  double E0 = pow(10.,par[2]);
  double NT = pow(10.,par[3]);

  double E   = var[0];
  double C  = pow((a-b)*E0/100.,a-b)*exp(b-a);
  if(E <= (a-b)*E0) return NT * pow(E/100.,a) * exp(-E/E0);
  return C* NT * pow(E/100.,b); // ph cm^(-2) s^(-1) keV^(-1)
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
  Nv->SetZTitle("N_{v} [ph/m^2/s/keV]");
  Nv->GetXaxis()->SetTitleOffset(1.5);
  Nv->GetYaxis()->SetTitleOffset(1.5);
  Nv->GetZaxis()->SetTitleOffset(1.2);
  Nv->GetXaxis()->CenterTitle();
  Nv->GetYaxis()->CenterTitle();
  Nv->GetZaxis()->CenterTitle();
  
  return Nv;
}



void MakeGRB()
{
}


void PlotGRB(double enph = 0,char name[100]="grb_65540.root")
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
  Ne->Scale(DT);
  Ne->SetYTitle("Ne [ph/keV/m^2]");
  Ne->SetXTitle(" Energy[keV] ");
  // Fv = Ne * de: [ph/s/m^2]
  TH1D *Fv = (TH1D*) Ne->Clone();
  Fv->SetTitle("Fv");
  Fv->SetName("Fv");
  Fv->SetYTitle("Fv [keV/keV/m^2]");
  
  // vFv = e * Ne * de: [keV/s/m^2]
  TH1D *vFv = (TH1D*) Ne->Clone();
  vFv->SetTitle("vFv");
  vFv->SetName("vFv");

  vFv->SetYTitle(" Flux ");
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
  SpectObj *sp = new SpectObj(Nv);              //ph/m² 
  
  std::pair<float,float> dir = std::make_pair((float)0.0,(float)0.0);
  
  sp->SaveParameters(0.0,dir);
  
  //////////////////////////////////////////////////
  TCanvas *clc = new TCanvas("clc","clc",600,800);
  TCanvas *csp = new TCanvas("csp","csp");
 
  csp->SetLogx();
  csp->SetLogy();
  
  clc->Divide(1,2);
  TH1D *Lct_GRB = sp->Integral_E(EMIN,EMAX);  // ph/m^2
  TH1D *Lct_GBM = sp->Integral_E(GBM1,GBM2);  // ph/m^2
 
  TH1D *Lct_LAT = sp->Integral_E(LAT1,LAT2);  // ph/m^2
  TH1D *Lct_EXT = sp->Integral_E(enph,EMAX);  // ph/m^2
  
  Lct_GRB->SetStats(0);
  Lct_LAT->SetStats(0);
  
  Lct_GBM->SetTitle("flux [ph]");
  Lct_GBM->SetXTitle("Time (s)");
  Lct_GBM->SetYTitle("photons/m^2");
    
  Lct_GRB->SetName("flux GRB [ph]");
  Lct_GBM->SetName("flux GBM [ph]");
  Lct_LAT->SetName("flux LAT [ph]"); 
  Lct_EXT->SetName("flux EXT [ph]"); 
 
  Lct_GRB->SetLineWidth(2);
  Lct_GBM->SetLineWidth(2);
  Lct_LAT->SetLineWidth(2);
  Lct_EXT->SetLineWidth(2);
  
  Lct_GRB->SetLineColor(2);
  Lct_GBM->SetLineColor(3);
  Lct_LAT->SetLineColor(4);
  Lct_EXT->SetLineColor(6);
 
  csp->cd();
  vFv->SetMinimum(.1);
  vFv->Draw("al");
  Fv->Draw("samel");
  Ne->Draw("samel");
  //  gDirectory->Delete("band");
  TF1 band("grb_f",Band,EMIN,1.0e+3,4); 
  //  band.SetParNames("a","b","Log10(E0)","Log10(Const)");
  band.SetLineStyle(2);
  band.SetParameters(-1.0,-2.5,1.0,4.0);
  band.SetParLimits(0,-2.0,0.0);
  band.SetParLimits(1,-4.0,-1.0);
  band.SetParLimits(2,1.0,3.0);
  band.SetParLimits(3,0.0,10.0);
  //  band.Draw("lsame");
  Ne->Fit("grb_f","v","lsame");
  std::cout<<"--------------------------------------------------"<<std::endl;
  std::cout<<" a     = "<<band.GetParameter(0)<<std::endl;
  std::cout<<" b     = "<<band.GetParameter(1)<<std::endl;
  std::cout<<" E0    = "<<pow(10.,band.GetParameter(2))<<std::endl;
  std::cout<<" Const = "<<pow(10.,band.GetParameter(3))<<std::endl;
  std::cout<<"--------------------------------------------------"<<std::endl;
  TLegend *leg = new TLegend(0.11,0.12,0.37,0.25);
  leg->SetFillStyle(0);
  leg->AddEntry(Ne," Ne  [ph/keV/m^2] ","lp");
  leg->AddEntry(Fv," Fv  [ph/m^2] ","lp");
  leg->AddEntry(vFv," vFv [keV/m^2] ","lp");
  leg->Draw();
  
  TH1D *Counts   = sp->GetSpectrum(); //ph
  Counts->SetName("Counts");
  TH1D *Lc = sp->GetTimes();  
  Lc->SetName("Counts[ph]");

  TF1  *f1 = new TF1("f1","[0]*1./x"); 
  TF1  *f2 = new TF1("f2","[0]*1./x^2");   
  
  
  if(ExtractPhotons)
    {  
      double time =0.0;
      double Interval;
      double energy;
      double flux;
      int i = 0;
      time = sp->interval(0.0,enph); // s
      while(time < TMAX && time>=0)
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
  double LowEnergyMax  = 1.05*Lct_GRB->GetMaximum();
  double HighEnergyMax = Lct_LAT->GetMaximum();
  double tmax = Lct_GRB->GetXaxis()->GetXmax();
  double tmin = Lct_GRB->GetXaxis()->GetXmin();
  /*
    Lct_LAT->Scale(LowEnergyMax/(1.1*HighEnergyMax));
    Lct_EXT->Scale(LowEnergyMax/(1.1*HighEnergyMax));
  */
  
  clc->cd(1);
  Lct_GRB->Draw("l");
  Lct_GBM->Draw("lsame");

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



  std::cout<<" Nph GRB ("<<EMIN<<","<<EMAX<<")  = "<<sp->Integral_T(Lct_GRB,0.0,TMAX)<<std::endl;
  std::cout<<" Nph GBM ("<<GBM1<<","<<GBM2<<")  = "<<sp->Integral_T(Lct_GBM,0.0,TMAX)<<std::endl;
  std::cout<<" Nph LAT ("<<LAT1<<","<<LAT2<<")  = "<<sp->Integral_T(Lct_LAT,0.0,TMAX)<<std::endl;
  if(ExtractPhotons) std::cout<<" Nph EXT ("<<enph<<","<<EMAX<<")  = "<<sp->Integral_T(Lct_EXT,0.0,TMAX)<<std::endl;
  std::cout<<" -------------------------------------------------- "<<std::endl;
  
  int iGBM1 = Ne->FindBin(GBM1);
  int iGBM2 = Ne->FindBin(GBM2);
  int iLAT1 = Ne->FindBin(LAT1);
  int iLAT2 = Ne->FindBin(LAT2);
  int iEXP = Ne->FindBin(enph);
  
  std::cout<<" Flux[ erg/cm^2] GRB ("<<EMIN<<","<<EMAX<<")  = "<<Fv->Integral(0,EBIN,"width")*1.0e-7/erg2meV<<" erg/cm^2"<<std::endl;
  std::cout<<" Flux[ erg/cm^2] GBM ("<<GBM1<<","<<GBM2<<")  = "<<Fv->Integral(iGBM1,iGBM2,"width")*1.0e-7/erg2meV<<" erg/cm^2"<<std::endl;
  std::cout<<" Flux[ erg/cm^2] LAT ("<<LAT1<<","<<LAT2<<")  = "<<Fv->Integral(iLAT1,iLAT2,"width")*1.0e-7/erg2meV<<" erg/cm^2"<<std::endl;
  if(ExtractPhotons) std::cout<<" Flux[ erg/cm^2] EXT ("<<enph<<","<<EMAX<<")  = "<<Fv->Integral(iEXP,EBIN,"width")*1.0e-7/erg2meV<<" erg/cm^2"<<std::endl;
  
}

int main(int argc, char** argv)
{
  
  std::cout<<" ****** GRB and ROOT test ****** "<<std::endl;
  std::string arg_name("");
  int current_arg = 1;
  long seed=-1;
  double enph=0.0;
  while(current_arg < argc)
    {
      arg_name = argv[current_arg];
      
      if("-extract"==arg_name)
	{
	  enph  = atof(argv[++current_arg]);
	}
      else if("-seed"==arg_name)
	{
	  seed=atoi(argv[++current_arg]);
	}
      current_arg++;
    }
  TApplication theApp("App",0,0);
  std::string paramFile = "$(GRBROOT)/src/test/GRBParam.txt";
  facilities::Util::expandEnvVar(&paramFile);
  
  params = new Parameters();  
  long seedN = params->GetGRBNumber();
  //  if(seed>=0) params->SetGRBNumber(seedN+seed);
  if(seed<1) seed=1;
  //////////////////////////////////////////////////
  params->ReadParametersFromFile(paramFile,seed);
  params->PrintParameters();
  GRBSim* m_grb = new GRBSim(params);
  //  m_grb->Fireball();
  m_grb->SaveNv((TH2D*)m_grb->Fireball());
  char name[100];
  sprintf(name,"grb_%d.root",params->GetGRBNumber());
  PlotGRB(enph,name);
  
  theApp.Run();
}
