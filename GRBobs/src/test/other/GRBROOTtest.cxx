#include <iostream>
//#include <vector>

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


//#include "../../GRBOBS/SpectObj.h"
#include "SpectObj/SpectObj.h"
#include "GRBobs/GRBobsSim.h"
#include "GRBobs/GRBobsConstants.h"

using namespace ObsCst;
double EMIN, EMAX, TMIN, TMAX, DT;
int    TBIN, EBIN;
//////////////////////////////////////////////////


double Band(double *var, double *par);
void help();
#define DEBUG=0

double Band(double *var, double *par)
{

  double a  = par[0];
  double b  = par[1];
  double E0 = pow(10.,par[2]);
  double NT = pow(10.,par[3]);

  double E   = var[0];
  double C   = pow((a-b)*E0,a-b)*exp(b-a);
  double H   = (a-b) * E0;
  if(E <= H) 
    return E * NT * pow(E,a) * exp(-E/E0);
  return E * C* NT * pow(E,b); // ph cm^(-2) s^(-1) keV^(-1)
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


void PlotGRB(double enph = 0,char name[100]="grb_65540.root")
{
  gStyle->SetCanvasColor(10);
  
  TCanvas *cNv = new TCanvas("cNv","Nv");
  cNv->SetLogy();
  cNv->SetLogz();
  std::cout<<" Loading "<<name<<std::endl;
  TH2D *Nv = Load(name); // [ph/keV/s/m²]
  
  // Ne =  [ph/keV/s/m²]
  gDirectory->Delete("Ne");
  Nv->ProjectionY("Ne");
  TH1D *Ne = (TH1D*) gDirectory->Get("Ne");

  Ne->Scale(DT);// Ne =  [ph/keV/m²]

  Ne->SetYTitle("Ne [ph/keV/m^{2}]");
  Ne->SetXTitle(" Energy[keV] ");
  // Fv = Ne * de: [ph/m^2]
  TH1D *Fv = (TH1D*) Ne->Clone();
  Fv->SetTitle("Fv");
  Fv->SetName("Fv");
  Fv->SetYTitle("Fv [keV/keV/m^{2}]");  
  // vFv = e * Ne * de: [keV/m^2]
  TH1D *vFv = (TH1D*) Ne->Clone();
  vFv->SetTitle("Fluxes");
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
  
  TCanvas *animatedCanvas = new TCanvas("animatedCanvas","animatedCanvas",500,500);
  animatedCanvas->SetLeftMargin(0.15);
  animatedCanvas->SetLogx();
  animatedCanvas->SetLogy();
  int i1=0;
  int i2;
  TH1D *AnimatedSpectrum;

  int l=1;
  for(int i=0; i < TBIN; i++)
    {
      i2=i;
      if(5==0)
	{
	  TString Name="frame";
	  AnimatedSpectrum =   Nv->ProjectionY("Animated",i1,i2);
	  AnimatedSpectrum->GetXaxis()->SetTitle("Energy [keV]");
	  AnimatedSpectrum->GetYaxis()->SetTitle("#nu F_{#nu} [erg/cm^{2}/s]");
	  AnimatedSpectrum->GetXaxis()->SetTitleOffset(1.1);
	  AnimatedSpectrum->GetYaxis()->SetTitleOffset(1.5);
	  AnimatedSpectrum->SetLineColor(2);
	  AnimatedSpectrum->SetLineWidth(2);
	  for(int j=0; j < EBIN; j++)
	    {
	      double ne = AnimatedSpectrum->GetBinContent(j+1);
	      double de = AnimatedSpectrum->GetBinWidth(j+1);
	      double ee = AnimatedSpectrum->GetBinCenter(j+1);
	      //	      AnimatedSpectrum->SetBinContent(i+1,de*ne);
	      AnimatedSpectrum->SetBinContent(j+1,ee*ne*de/624151.0*1e-7);
	    }  
	  AnimatedSpectrum->SetMaximum(1e-3);
	  AnimatedSpectrum->SetMinimum(1e-10);
	  AnimatedSpectrum->Draw("l");
	  animatedCanvas->Update();
	  Name+=l;
	  Name+=".gif";
	  animatedCanvas->Print(Name);
	  l++;

	  i1=i;
	}
    }
  
  

  //////////////////////////////////////////////////
  cNv->cd();
  Nv->Draw("surf");
  
  bool ExtractPhotons = true;
  if(enph<=0) ExtractPhotons = false;
  else if (enph<EMIN)  enph = EMIN;
  
  //////////////////////////////////////////////////
  SpectObj *sp = new SpectObj(Nv,0);              //ph

  std::pair<float,float> dir = std::make_pair((float)0.0,(float)0.0);
  
  //  sp->SaveParameters(0.0,dir);
  
  //////////////////////////////////////////////////
  TCanvas *clc = new TCanvas("clc","clc",600,800);
  TCanvas *csp = new TCanvas("csp","csp");
 
  csp->SetLogx();
  csp->SetLogy();
  
  clc->Divide(1,2);
  TH1D *Lct_GBM = sp->Integral_E(GBM1,GBM2);  // ph
  TH1D *Lct_BATSE1 = sp->Integral_E(BATSE1,BATSE2);  // ph
  TH1D *Lct_BATSE2 = sp->Integral_E(BATSE2,BATSE3);  // ph
  TH1D *Lct_BATSE3 = sp->Integral_E(BATSE3,BATSE4);  // ph
  TH1D *Lct_BATSE4 = sp->Integral_E(BATSE4,BATSE5);  // ph
 
  TH1D *Lct_LAT = sp->Integral_E(LAT1,LAT2);  // ph
  TH1D *Lct_EXT = sp->Integral_E(enph,EMAX);  // ph
  
  Lct_GBM->SetStats(0);
  Lct_LAT->SetStats(0);
  
  Lct_GBM->SetNameTitle("flux [ph]","flux [ph]");
  Lct_GBM->SetXTitle("Time (s)");
  Lct_GBM->SetYTitle("photons");
  
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
  Lct_GBM->SetLineWidth(1);
  Lct_LAT->SetLineWidth(2);
  Lct_EXT->SetLineWidth(1);
  Lct_EXT->SetLineStyle(2);
  
  Lct_BATSE1->SetLineColor(2);
  Lct_BATSE2->SetLineColor(3);
  Lct_BATSE3->SetLineColor(4);
  Lct_BATSE4->SetLineColor(6);
  //  Lct_GBM->SetLineColor(1);
  Lct_LAT->SetLineColor(4);
  //  Lct_EXT->SetLineColor(6);
 
  csp->cd();
  vFv->SetMinimum(.1);
  vFv->Draw("al");
  Fv->Draw("samel");
  Ne->Draw("samel");
  //  gDirectory->Delete("band");
  
  TF1 band("grb_f",Band,EMIN,1.0e+4,4); 
  TF1 pl("PL_f","10^([0])*(x/10000)^(-[1])",3.0e+3,1.0e+6);
 
  band.SetParNames("a","b","Log10(E0)","Log10(Const)");
  band.SetLineStyle(2);
  pl.SetLineStyle(2);
  pl.SetLineColor(2);

  band.SetParLimits(0,-2.0,1.0);
  band.SetParLimits(1,-3.0,-1.0);
  band.SetParLimits(2,1.0,3.0);
  band.SetParLimits(3,0.0,10.0);

  band.SetParameter(0,-0.86);
  //  band.FixParameter(1,-2.125);
  band.SetParameter(2,2.5);
  band.SetParameter(3,3.0);
  //band.Draw("lsame");
  pl.SetParameters(4.0,-2.5);
  pl.SetParLimits(0,-8,8);
    
  Fv->Fit("grb_f","QR+","lsame");

  double a=band.GetParameter(0);
  double b=band.GetParameter(1);
  double E0=pow(10.,band.GetParameter(2));
  double Ep=(a+2)*E0;
  std::cout<<"-- BAND FUNCTION Fv FIT ---------------------------"<<std::endl;
  std::cout<<" a     = "<<a<<std::endl;
  std::cout<<" b     = "<<b<<std::endl;
  std::cout<<" E0    = "<<E0<<std::endl;
  std::cout<<" Ep    = "<<Ep<<std::endl;
  std::cout<<" Const = "<<pow(10.,band.GetParameter(3))<<std::endl;
  std::cout<<"--------------------------------------------------"<<std::endl;
  //////////////////////////////////////////////////
  Ne->Fit("PL_f","QR+","lsame");
  std::cout<<"---PL FUNTION Nv FIT------------------------------"<<std::endl;
  std::cout<<" No    = "<<pow(10.,pl.GetParameter(0))<<std::endl;
  std::cout<<" Index = "<<pl.GetParameter(1)<<std::endl;
  std::cout<<"--------------------------------------------------"<<std::endl;

  Fv->Fit("PL_f","QR+","lsame");
  std::cout<<"---PL FUNTION Fv FIT------------------------------"<<std::endl;
  std::cout<<" No    = "<<pow(10.,pl.GetParameter(0))<<std::endl;
  std::cout<<" Index = "<<pl.GetParameter(1)<<std::endl;
  std::cout<<"--------------------------------------------------"<<std::endl;

  vFv->Fit("PL_f","QR+","lsame");
  std::cout<<"---PL FUNTION vFv FIT------------------------------"<<std::endl;
  std::cout<<" No    = "<<pow(10.,pl.GetParameter(0))<<std::endl;
  std::cout<<" Index = "<<pl.GetParameter(1)<<std::endl;
  std::cout<<"--------------------------------------------------"<<std::endl;
   
  TLegend *leg = new TLegend(0.11,0.12,0.37,0.25);
  leg->SetFillStyle(0);
  leg->AddEntry(Ne," Ne  [ph/keV/m^{2}] ","lp");
  leg->AddEntry(Fv," F_{#nu}  [ph/m^{2}] ","lp");
  leg->AddEntry(vFv," #nu F_{#nu} [keV/m^{2}] ","lp");
  leg->Draw();
  
  TH1D *Counts   = sp->CloneSpectrum(); //ph
  TH1D *Lc = sp->CloneTimes();  
  Lc->SetTitle("Counts[ph]");
  Lc->SetXTitle("Time (s)");
  Lc->SetYTitle("photons");

  TF1  *f1 = new TF1("f1","[0]*1./x"); 
  TF1  *f2 = new TF1("f2","[0]*1./x^2");   
  
  
  if(ExtractPhotons)
    {  
      double time = 0.0;
      double Interval;
      double energy;
      double flux;
      int i = 0; 
      Interval = sp->interval(time,enph); // s
      time+=Interval;
      while(time < TMAX && time>=0)
	{
	  energy   = sp->energy(time,enph);   // keV
	  flux     = sp->flux(time,enph);     // ph/s/m2
	  Counts->Fill(energy); // ph
	  Lc->Fill(time);       // ph
	  Interval = sp->interval(time,enph); // s  
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
      Counts->Scale(1./6.0);
      Counts->Draw("E1same");
      Counts->SetStats(1);
    }
  //////////////////////////////////////////////////      
  //clc->cd();
  double LowEnergyMax  = 1.05*Lct_GBM->GetMaximum();
  double HighEnergyMax = Lct_LAT->GetMaximum();
  double tmax = Lct_GBM->GetXaxis()->GetXmax();
  double tmin = Lct_GBM->GetXaxis()->GetXmin();
  /*
    Lct_LAT->Scale(LowEnergyMax/(1.1*HighEnergyMax));
    Lct_EXT->Scale(LowEnergyMax/(1.1*HighEnergyMax));
  */
  
  clc->cd(1);
  Lct_GBM->Draw("l");
  Lct_BATSE1->Draw("lsame");
  Lct_BATSE2->Draw("lsame");
  Lct_BATSE3->Draw("lsame");
  Lct_BATSE4->Draw("lsame");
  TLegend *legL = new TLegend(0.7,0.7,0.99,0.99);
  legL->AddEntry(Lct_BATSE1,"BATSE (ch1)");
  legL->AddEntry(Lct_BATSE2,"BATSE (ch2)");
  legL->AddEntry(Lct_BATSE3,"BATSE (ch3)");
  legL->AddEntry(Lct_BATSE4,"BATSE (ch4)");
  legL->AddEntry(Lct_GBM,"GBM (tot)");
  legL->Draw();
  clc->cd(2);
  TLegend *legH = new TLegend(0.7,0.8,0.99,0.99);
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
      std::cout<<" flux (exp)        : = "<<Counts->Integral(0,EBIN,"width")*1.0e-7/erg2meV/6.0<<" erg/cm^2"<<std::endl;
      
    } 
  else
    {
      Lct_LAT->Draw("l");
    }

  legH->AddEntry(Lct_LAT,"LAT (tot)");
  legH->Draw();


  std::cout<<" Nph BATSE1 ("<<BATSE1<<","<<BATSE2<<")  = "<<sp->Integral_T(Lct_BATSE1,0.0,TMAX)<<std::endl;
  std::cout<<" Nph BATSE1 ("<<BATSE2<<","<<BATSE3<<")  = "<<sp->Integral_T(Lct_BATSE2,0.0,TMAX)<<std::endl;
  std::cout<<" Nph BATSE1 ("<<BATSE3<<","<<BATSE4<<")  = "<<sp->Integral_T(Lct_BATSE3,0.0,TMAX)<<std::endl;
  std::cout<<" Nph BATSE1 ("<<BATSE4<<","<<BATSE5<<")  = "<<sp->Integral_T(Lct_BATSE4,0.0,TMAX)<<std::endl;
  std::cout<<" Nph GBM ("<<GBM1<<","<<GBM2<<")  = "<<sp->Integral_T(Lct_GBM,0.0,TMAX)<<std::endl;
  std::cout<<" Nph LAT ("<<LAT1<<","<<LAT2<<")  = "<<sp->Integral_T(Lct_LAT,0.0,TMAX)<<std::endl;
  if(ExtractPhotons) std::cout<<" Nph EXT ("<<enph<<","<<EMAX<<")  = "<<sp->Integral_T(Lct_EXT,0.0,TMAX)<<std::endl;
  std::cout<<" -------------------------------------------------- "<<std::endl;
  
  /*int iBATSE1 = Ne->FindBin(BATSE1);
  int iBATSE2 = Ne->FindBin(BATSE2);
  int iGBM1 = Ne->FindBin(GBM1);
  int iGBM2 = Ne->FindBin(GBM2);
  int iLAT1 = Ne->FindBin(LAT1);
  int iLAT2 = Ne->FindBin(LAT2);
  int iEXP = Ne->FindBin(enph);
  */
  double fBATSE1 = sp->GetFluence(BATSE1,BATSE2);
  double fBATSE2 = sp->GetFluence(BATSE2,BATSE3);
  double fBATSE3 = sp->GetFluence(BATSE3,BATSE4);
  double fBATSE4 = sp->GetFluence(BATSE4,BATSE5);
  double fLAT   = sp->GetFluence(LAT1,LAT2);
  double fGBM   = sp->GetFluence(GBM1,GBM2);
  double fEXP   = sp->GetFluence(enph,emax);
  double fTOT   = sp->GetFluence();
  Ep = vFv->GetBinCenter(vFv->GetMaximumBin());
  std::cout<<"* Theoretical values:  *****************************"<<std::endl;
  std::cout<<" T90 = "<<sp->GetT90()<<" Epeak = "<<Ep<<std::endl;
  std::cout<<" GRB   flux ("<< emin <<","<< emax <<") = "<<fTOT<<" erg/cm^2"<<std::endl;
  std::cout<<" BASTE flux (ch1) ("<<BATSE1<<","<<BATSE2<<") = "<<fBATSE1<<" erg/cm^2"<<std::endl;
  std::cout<<" BASTE flux (ch2) ("<<BATSE2<<","<<BATSE3<<") = "<<fBATSE2<<" erg/cm^2"<<std::endl;
  std::cout<<" BASTE flux (ch3) ("<<BATSE3<<","<<BATSE4<<") = "<<fBATSE3<<" erg/cm^2"<<std::endl;
  std::cout<<" BASTE flux (ch4) ("<<BATSE4<<","<<BATSE5<<") = "<<fBATSE4<<" erg/cm^2"<<std::endl;
  std::cout<<" GBM   flux ("<< GBM1 <<","<< GBM2 <<") = "<<fGBM<<" erg/cm^2"<<std::endl;
  std::cout<<" LAT   flux ("<< LAT1 <<","<< LAT2 <<") = "<<fLAT<<" erg/cm^2"<<std::endl;
  if(ExtractPhotons) std::cout<<" EXP   flux ("<< enph <<","<< emax <<") = "<<fEXP<<" erg/cm^2"<<std::endl;
  
  /*
  std::cout<<" Flux[ erg/cm^2] GRB ("<<EMIN<<","<<EMAX<<")  = "<<Fv->Integral(0,EBIN,"width")*1.0e-7/erg2meV<<" erg/cm^2"<<std::endl;
std::cout<<" Flux[ erg/cm^2] BATSE ("<<BATSE1<<","<<BATSE2<<")  = "<<Fv->Integral(iBATSE1,iBATSE2,"width")*1.0e-7/erg2meV<<" erg/cm^2"<<std::endl;
  std::cout<<" Flux[ erg/cm^2] GBM ("<<GBM1<<","<<GBM2<<")  = "<<Fv->Integral(iGBM1,iGBM2,"width")*1.0e-7/erg2meV<<" erg/cm^2"<<std::endl;
  std::cout<<" Flux[ erg/cm^2] LAT ("<<LAT1<<","<<LAT2<<")  = "<<Fv->Integral(iLAT1,iLAT2,"width")*1.0e-7/erg2meV<<" erg/cm^2"<<std::endl;
  if(ExtractPhotons) std::cout<<" Flux[ erg/cm^2] EXT ("<<enph<<","<<EMAX<<")  = "<<Fv->Integral(iEXP,EBIN,"width")*1.0e-7/erg2meV<<" erg/cm^2"<<std::endl;
  */
}


void MakeGRB(int NGRB=1, double enph=0, bool gbm=false)
{
  std::cout<<" ****** GRB and ROOT test ****** "<<std::endl;
  
  std::string path = ::getenv("GRBOBSROOT");
  std::string paramFile = path+"/src/test/GRBParam.txt";
  
  GRBobsParameters *params = new GRBobsParameters();  
  //////////////////////////////////////////////////
  params->ReadParametersFromFile(paramFile,NGRB);
  GRBobsSim* m_grb = new GRBobsSim(params);
  m_grb->MakeGRB();
  m_grb->SaveNv();
  if (gbm)  m_grb->GetGBMFlux();

  char name[100];
  sprintf(name,"grbobs_%d.root",params->GetGRBNumber());
  delete m_grb;
  //  delete params; ??
  PlotGRB(enph,name);
  
 }



int main(int argc, char** argv)
{
  
  std::string arg_name("");
  int current_arg = 1;
  double enph=0.0;
  int ngrb=1;
  bool video_out=true;
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
      current_arg++;
    }
  
  TApplication theApp("App",0,0);

  MakeGRB(ngrb,enph,gbm);
  theApp.Run();
}

void help()
{
  std::cout<<"--------------------------------------------------"<<std::endl;
  std::cout<<"   GRB ROOT test program author: Nicola Omodei    "<<std::endl;
  std::cout<<"   Option are: "<<std::endl;
  std::cout<<"  -extract [enph] : etxtract photons above enph"<<std::endl;
  std::cout<<"  -grb [N] processes the N grb in the $GRBROOT/src/test/GRBParam.txt file "<<std::endl;
  std::cout<<"  -gbm      Fit the GBM spectrum with as Band function, as afunction of the time"<<std::endl;
  std::cout<<"--------------------------------------------------"<<std::endl;

}

