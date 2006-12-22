/*!\file GRBTest.cxx
 * \brief test program for GRB simulation studies.
 * 
 * This executable uses ROOT to display several histograms at run time.
 * To use it, type on the prompt ./GRBTest.exe
 * After some output, an empty canvas named \b c1 will appear, 
 * and the program hangs at the line:
 * "Type a time (in sec) or 0 for the complete evolution"
 * \arg A time \c t entered will provoke the computation of the flux at \c t
 * seconds after initiation of the burst in the GLAST referential frame.
 * Then, the canvas \c1 will show on top the number of photons reaching the 
 * detector as a function of their energy (in log scale eV), and on bottom the
 * sampled photons from the preceding curve (note that the lowest energy for 
 * sampling is set to 10 Mev.)
 * \arg If 0 is typed, the complete simulation of GRB photon arrival to the 
 * detector is simulated. The canvas displays the continuous evolution of 
 * the spectrum. At the end, 4 new canvas pops up, showing several summary 
 * histograms.
 */
#include <iterator>
#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>
//Include files for spectrum...

//include files for ROOT...

#include "TROOT.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom.h"
#include "TTree.h"
#include "TBranch.h"
#include "TObjArray.h"
//#include "BranchObject.h"
#include "TFile.h"

//GRB include files...
#include "../../src/GRB/GRBShell.h"
#include "../../src/GRB/GRBShock.h"
#include "../../src/GRB/GRBConstants.h"
#include "../../src/GRB/GRBSim.h"
//#include "../../src/GRB/GRBSpectrum.h"
#include "CLHEP/Random/RandFlat.h"
using namespace std;

//! The size of the energy array is taken from GRBConstant.h file in FluxSvc
#define ENERGYSIZE cst::enstep
//! The size of the time array is taken from GRBConstant.h file in FluxSvc
#define TIMESIZE cst::nstep
/*! the flux is a bidimentional array (matrix?) funtion of energy and time.
 * Its Unities are 
 *
 */
const double TIME=10.0;
const double EVENTS=100000;
double m_flux[ENERGYSIZE][TIMESIZE];
double timex[TIMESIZE];
double energyx[ENERGYSIZE];

static const char * default_arg="GRBSpectrum";
//int test1(int argc, char** argv);
void test2();
//void Generate();

void help() {
   std::cout << 
      "   Simple test program for Transient Sources.\n"
      "   Command line args are \n"
      "      '-events <number of events to create>'\n"
      "      '-time <time in seconds>'    for the maximum time\n"
      "      '-list' lists the available spectra\n"
      "      '-help' for this help"
      << std::endl;
}

int main()
{
  //  test1(argc,argv);
  test2();
  return 0;
}

double CalculateFluence(double ee,double e1=cst::enmin*1e-9,double e2=cst::enmax*1e-9)
{
  double fluence;
  if (ee>=e1 && ee<=e2) 
    {
      fluence=ee;
    }
  else 
    {
      fluence=0.0;
    }
  return fluence;
}

void test2()
{
  int flag=1;
  std::vector<double> spectrum,energy,deltae;
  
  GRBSim* _myGRB = new GRBSim();
  /// This initializes the GRB simulation...
  cout<<"******     Initializing ROOT    ******"<<endl;
  static const   TROOT root("GRB", "GRB Simulator");
  _myGRB->Start();
  
  double tmax=_myGRB->Tmax();
  double m_dt=tmax/cst::nstep;
  
  energy=_myGRB->Energy();
  deltae=_myGRB->DeltaE();
  
  TApplication theApp("App", 0, 0);
  
  float phenergy;
  
  TH1D *hh1 = new TH1D("hh1","hh1",cst::enstep,log10(cst::enmin),log10(cst::enmax));
  TH1D *hh2 = new TH1D("hh2","hh2",cst::enstep,log10(cst::enmin),log10(cst::enmax));
  
  double ph;
  double eextra;
  double tim;
  int temp;
  double etot2;
  temp=1;
  
  TCanvas *cc1 = new TCanvas();  
  cc1->Divide(0,2);
  cout<<"Type a time (in sec) or a 0 for the complete evolution:  "<<endl;
  while(temp==1)
    {
      cin>>tim;
      hh1->Reset();
      hh2->Reset();
      spectrum.clear();
      
      if (tim>0)
	{
	  etot2=0.0;
	  _myGRB->ComputeFlux(tim);
	  spectrum=_myGRB->Spectrum();
	  
	  for (int en=0;en<cst::enstep;en++)
	    {
	      if (energy[en]<cst::enph) 
		{ 
		  hh1->SetBinContent(en+1,1.0);
		}
	      else
		{
		  hh1->SetBinContent(en+1,spectrum[en]);
		}
	    }

	  cout<<"Time  = "<< tim <<" Ftot  [eV/s/m^2]= "<<_myGRB->IFlux(cst::enph)<<endl;
	  cout<<"                    Counts[ph/s/m^2]=" <<_myGRB->IRate(cst::enph)<<endl;
	  cout<<"                    Energy[eV/m^2]=" <<_myGRB->IEnergy(cst::enph)<<endl;
	  
	  cc1->cd(1);
	  hh1->Draw("AL");
	  cc1->Update(); 
    
	  cc1->cd(2);
	  hh1->Draw("AL");
	  cc1->Update();     
	  //      hh1->Integral(e1);
	  eextra=0.0;
	  int j=0;
	  while (eextra<_myGRB->IEnergy(cst::enph))
	    {
	      phenergy=_myGRB->DrawPhotonFromSpectrum(spectrum,RandFlat::shoot(1.0),cst::enph);
	      ph=log10(phenergy*1e+9);
	      hh2->Fill(ph);
	      eextra+=phenergy*1.0e+9;
	      j++;
	    }
	  
	  cout<< "Number of photons extracted: "<<j<<" Energy extracted = "<<eextra<<endl<<endl;
	  cc1->cd(2);
	  hh2->Draw();
	  cc1->Update();
	} else 
	  {
	    tim=-1.0;
	    temp=0;
	  }
    }
  cc1->Clear();
  cc1->Delete();
  cout<<temp<<endl;

  ////////////////  Compute the Total Flux ////////////////
  //------------------------------------------------------//  
  double ftottot=0.0;
  std::vector<double> m_time;
  TH1D *h1 = new TH1D("h1","h1",cst::enstep,log10(cst::enmin),log10(cst::enmax));
  
  TCanvas *cc2 = new TCanvas();
  h1->GetXaxis()->SetTitle("Log(Energy [eV])");
  h1->GetXaxis()->SetLabelSize(0.05);
  h1->GetXaxis()->SetTitleOffset(1.5);
  //  h1->GetYaxis()->SetTitle("flux[ph/s/eV/m^2]");  
  h1->GetYaxis()->SetTitle("flux[erg/s/cm^2]");  
  //  h1->SetFillColor(5);
  
  for (int t=0;t<cst::nstep;t++)
    {
      for (int en=0;en<cst::enstep;en++)
	{
	  m_flux[en][t]=0.0;
	  cout<<energy[en]<<"  "<<deltae[en]
	      <<"e/de"<<energy[en]/deltae[en]<<endl;
	}
    }
  for (int t=0;t<cst::nstep;t++)
    {
      m_time.push_back(t*m_dt+m_dt);
      // Compute the flux @ time
      spectrum.clear();
      _myGRB->ComputeFlux(m_time[t]);    
      spectrum=_myGRB->Spectrum(); // is in photons/s/MeV/m^2
      for (int en=0;en<cst::enstep;en++)
	{
	  m_flux[en][t]=spectrum[en]*energy[en]/(cst::erg2MeV*1.0e+6);  
	  // m_flux is in erg/s/MeV/m^2
	}
      
      if (flag==1) 
	{
	  for (int en=0;en<cst::enstep;en++)
	    {
	      h1->SetBinContent(en+1,m_flux[en][t]*1.0e-4); //erg/cm^2/s/MeV
	      //	h1->SetBinContent(en+1,Flux(en)); //ph/m^2/s/eV
	    }
	  cc2->cd();
	  h1->Draw("AL");
	  h1->GetXaxis()->SetTitle("Log(Energy[eV])");
	  h1->GetXaxis()->SetTitleSize(0.035);
	  h1->GetXaxis()->SetLabelSize(0.035);
	  h1->GetXaxis()->SetTitleOffset(1.4);
	  h1->GetYaxis()->SetTitle("Fv[erg/cm^2/s/MeV]");
	  h1->GetYaxis()->SetTitleSize(0.035);
	  h1->GetYaxis()->SetLabelSize(0.035);
	  h1->GetYaxis()->SetTitleOffset(1.4);   
	  
	  h1->SetLineWidth(3);
	  h1->SetLineColor(4);
	  cc2->Update();   
	}
      ftottot+=_myGRB->IEnergy(cst::enmin)/(cst::erg2MeV*1.0e+10)*_myGRB->Area();
      //erg
      cout<<"Time / tmax = "<<m_time[t]<<"/" <<tmax<<" Ftot [erg]= "<<ftottot<<endl;
    }
  cout<<"Fluence [erg/cm^2] = "<<ftottot/_myGRB->Area()<<endl;
  
  //////////////////////////////////////////////////////////
  for (int t=0;t<cst::nstep;t++)
    {
      timex[t]=m_time[t];
    }
  for (int en=0;en<=cst::enstep;en++)
    {
      energyx[en]=energy[en];
    }
  
  TCanvas *c1 = new TCanvas("c1","GRB Flux");
  TCanvas *c1b = new TCanvas("c1b","Photons Count");
  TCanvas *c2 = new TCanvas("c2","GRB Light Curve");
  TCanvas *c4 = new TCanvas("c4","Countour Plot");
  TCanvas *c5 = new TCanvas("c5","3D Flux Rapresentation");
  //cout<<"xmin = "<< energyx[0] << "   xmax = "<<energyx[cst::nstep-1]<< "  ymin = "<<timex[0] <<"   ymax =" << timex[cst::nstep-1] <<endl;
  
  c1->SetLogx();
  c1->SetLogy();
  c1->SetGridx();
  c1->SetGridy();
  
  c1b->SetLogx();
  c1b->SetLogy();
  c1b->SetGridx();
  c1b->SetGridy();
  
  
  TGraph *gfl = new TGraph(cst::enstep);
  TGraph *gph = new TGraph(cst::enstep);
  
  //TGraph *glc = new TGraph(cst::nstep);
  TGraph *gch1 = new TGraph(cst::nstep);
  TGraph *gch2 = new TGraph(cst::nstep);
  TGraph *gch3 = new TGraph(cst::nstep);
  TGraph *gch4 = new TGraph(cst::nstep);
  
  TH2D *h2 = new TH2D("h2","H2",cst::enstep-2,log10(energyx[1]),log10(energyx[cst::enstep]),cst::nstep-1,timex[0],timex[cst::nstep-1]);
  //TH2D *h2 = new TH2D("h2","H2",cst::enstep-1,2,20,cst::nstep-1,0,1);
  //h1->SetOption("B");
  
  double Fv[cst::enstep];
  double ph_m2[cst::enstep];
  double lct[cst::nstep];
  double ch1[cst::nstep];
  double ch2[cst::nstep];
  double ch3[cst::nstep];
  double ch4[cst::nstep];
  
  for (int t=0;t<cst::nstep;t++)
    {
      ch1[t]=0.0;
      ch2[t]=0.0;
      ch3[t]=0.0;
      ch4[t]=0.0;
      lct[t]=0.0;
    }
  double loge[cst::enstep];
  for (int en=0;en<cst::enstep;en++)
    {
      Fv[en]=0.0;
      ph_m2[en]=0.0;
      loge[en]=log10(energyx[en]);    
      for (int t=0;t<cst::nstep;t++) 
	{
	  m_flux[en][t]<=1e-25?(m_flux[en][t]=1.0e-25):(m_flux[en][t]);
	  // m_flux is in erg/s/MeV/m^2
	  Fv[en]+=(m_flux[en][t])*m_dt;  
	  /// is in erg/m^2/MeV
	  ph_m2[en]+=(m_flux[en][t])*(cst::erg2MeV*1.0e+6)/
	    energyx[en]*m_dt;
	  // is in ph/MeV/m^2
	  lct[t]+=(m_flux[en][t])*(cst::erg2MeV)
	    /energyx[en]
	    *m_dt; /// is in ph/m^2
      
	  if (energyx[en]<cst::ch1H &&energyx[en]>=cst::ch1L)
	    {
	      ch1[t]+=(m_flux[en][t])*(cst::erg2MeV)
		/energyx[en]
		*m_dt; /// is in ph/m^2
	    }
	  if (energyx[en]<cst::ch2H &&energyx[en]>=cst::ch2L)
	    {
	      ch2[t]+=(m_flux[en][t])*(cst::erg2MeV)
		/energyx[en]
		*m_dt; /// is in ph/m^2
	    }
	  if (energyx[en]<cst::ch3H &&energyx[en]>=cst::ch3L)
	    {
	      ch3[t]+=(m_flux[en][t])*(cst::erg2MeV)
		/energyx[en]
		*m_dt; /// is in ph/m^2
	    }
	  if (energyx[en]<cst::ch4H &&energyx[en]>=cst::ch4L)
	    {
	      ch4[t]+=(m_flux[en][t])*(cst::erg2MeV)
		/energyx[en]
		*m_dt; /// is in ph/m^2
	    }
	  h2->Fill(loge[en],timex[t],m_flux[en][t]);
	  //      cout<<log10(energyx[en])<<timex[t]<<log10(m_flux[en][t])<<endl;
	}
    }
  
  gfl=new TGraph(cst::enstep,energyx,Fv);
  gph=new TGraph(cst::enstep,energyx,ph_m2);
  
  //glc=new TGraph(cst::nstep,timex,lct);
  
  gch1=new TGraph(cst::nstep,timex,ch1);
  gch2=new TGraph(cst::nstep,timex,ch2);
  gch3=new TGraph(cst::nstep,timex,ch3);
  gch4=new TGraph(cst::nstep,timex,ch4);
  
  gfl->SetLineWidth(2);
  //glc->SetLineWidth(2);
  gfl->SetLineColor(4);
  
  gch1->SetLineColor(2);
  gch2->SetLineColor(3);
  gch3->SetLineColor(4);
  gch4->SetLineColor(5);
  
  
  c1->cd();
  gfl->Draw("ALP");
  gfl->GetXaxis()->SetTitle("Energy [eV]");
  gfl->GetXaxis()->SetTitleSize(0.035);
  gfl->GetXaxis()->SetLabelSize(0.035);
  gfl->GetXaxis()->SetTitleOffset(1.4);
  gfl->GetYaxis()->SetTitle("Fv[erg/cm^2/MeV]");
  gfl->GetYaxis()->SetTitleSize(0.035);
  gfl->GetYaxis()->SetLabelSize(0.035);
  gfl->GetYaxis()->SetTitleOffset(1.4);
  gfl->Draw("ALP");
  
  c1b->cd();
  gph->Draw("ALP");
  gph->GetXaxis()->SetTitle("Log(Energy [eV])");
  gph->GetXaxis()->SetTitleSize(0.035);
  gph->GetXaxis()->SetLabelSize(0.035);
  gph->GetXaxis()->SetTitleOffset(1.4);
  gph->GetYaxis()->SetTitle("Particle/m^2/MeV]");
  gph->GetYaxis()->SetTitleSize(0.035);
  gph->GetYaxis()->SetLabelSize(0.035);
  gph->GetYaxis()->SetTitleOffset(1.4);
  gph->Draw("ALP");
  
  c1b->Update();
  
  char legendName1[100];
  char legendName2[100];
  char legendName3[100];
  char legendName4[100];
  //char legendName5[100];
  
  c2->cd();
  gch1->Draw("ALP");
  gch1->SetTitle("Light Curves");
  gch1->GetXaxis()->SetTitle("Time [sec]");
  gch1->GetXaxis()->SetTitleSize(0.035);
  gch1->GetXaxis()->SetLabelSize(0.035);
  gch1->GetXaxis()->SetTitleOffset(1.4);
  gch1->GetYaxis()->SetTitle("Counts [ph/m^2]");
  gch1->GetYaxis()->SetTitleSize(0.035);
  gch1->GetYaxis()->SetLabelSize(0.035);
  gch1->GetYaxis()->SetTitleOffset(1.4);
  //glc->Draw("ALP");
  gch1->Draw("ALP");
  gch2->Draw("LP");
  gch3->Draw("LP");
  gch4->Draw("LP");

  
  TLegend *leg = new TLegend(0.6,0.75,0.98,0.98);
  leg->SetTextFont(11);
  
  sprintf(legendName1,"(%6g GeV - %6g GeV) ",cst::ch1L*1.0e-9,cst::ch1H*1.0e-9);
  sprintf(legendName2,"(%6g GeV - %6g GeV) ",cst::ch2L*1.0e-9,cst::ch2H*1.0e-9);
  sprintf(legendName3,"(%6g GeV - %6g GeV) ",cst::ch3L*1.0e-9,cst::ch3H*1.0e-9);
  sprintf(legendName4,"(%6g GeV - %6g GeV) ",cst::ch4L*1.0e-9,cst::ch4H*1.0e-9);
  //  sprintf(legendName5,"(%6g GeV - %6g GeV) ",cst::ch1L,cst::ch1H);

  leg->AddEntry(gch1,legendName1,"l");
  leg->AddEntry(gch2,legendName2,"l");
  leg->AddEntry(gch3,legendName3,"l");
  leg->AddEntry(gch4,legendName4,"l");
  //  leg->AddEntry(gch5,legendName5,"l");
  //leg->AddEntry(glc,"Sum","l");
  leg->Draw();
  
  
  c2->Update();
  
  c4->cd(); 
  h2->Draw("CONT");
  c4->Update();
  
  c5->cd();
  h2->Draw("surf");
  h2->GetXaxis()->SetTitle("Log(Energy [eV])");
  //  h2->GetYaxis()->SetTitleSize(0.05);
  //  h2->GetXaxis()->SetLabelSize(0.05);
  h2->GetXaxis()->SetTitleOffset(2);
  h2->GetYaxis()->SetTitle("Time[sec]");
  //h2->GetYaxis()->SetTitleSize(0.035);
  //  h2->GetYaxis()->SetLabelSize(0.05);
  h2->GetYaxis()->SetTitleOffset(2);
  h2->GetZaxis()->SetTitle("vFv[erg/s/m^2/MeV]");
// h2->GetYaxis()->SetTitleSize(0.035);
h2->GetZaxis()->SetLabelSize(0.05);
  h2->GetZaxis()->SetTitleOffset(1.3);
  h2->Draw("surf");
  c5->Update();
  ////////////////////////////////////////
  delete _myGRB;
  theApp.Run();
  cin>>tim;
}

