// Description:
//    GRBROOTTest is a test program for GRB simulation studies.
//    This executable uses ROOT to display several histograms at run time.
//    It is just an interface to GRBSim, which actually does all the job.
//
// Usage:
//    test_GRBROOT.exe [options] [option_args]
//
// Options:
//    -help or -h print all the options
//    -time <time in sec> fix the maximum time
//    -save Save the parameters for the model and the output in a file ("GRBdata.txt").
//    -root Save the histograms in a in a root file ("histos.root")
//    -seed <seed> set the seed for the random engine
//    -mute It turns off the graphical output
//    -nomovie It turn off the visualization of the evolution of the flux
//    -cycle <N> it runs the program N times
//    -gif  It save the moovie of the evolving flux in a gif file.
//
// Graphical output:
//    After a period of initialization, a canvas pops up showing the 
//    complete evolution of the burst in function of the time. 
//    At the end, 4 new canvas pops up, showing several summary 
//    histograms.
//
// Authors
//    Nicola Omodei
//    Johann Cohen-Tanugi
//
// [Header File]
#include <iterator>
#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>
// include facilities
#include "facilities/Util.h"
//include files for ROOT...
#include "TROOT.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom.h"
#include "TTree.h"
#include "TBranch.h"
#include "TObjArray.h"
#include "TFile.h" 
//GRB include files...
#include "../../GRB/GRBShell.h"
#include "../../GRB/GRBShock.h"
#include "../../GRB/GRBConstants.h"
#include "../../GRB/GRBSim.h"
// CLHEP include
#include "CLHEP/Random/RandFlat.h"
using namespace std;

// The size of the energy array is taken from GRBConstant.h
#define ENERGYSIZE cst::enstep
// The size of the time array is taken from GRBConstant.h
#define TIMESIZE cst::nstep
////////////////////////////////////////////////////////////
// Channel for the output Light Curves, all the values are in eV
// the values are: ch*L : channel * lowest energy (eV)
//                 ch*H : channel * highest enery (eV)
//                 Nch* : channel * noise value (counts)
// GBM 1 Channel:
const double ch1L      =   50.0e+3;
const double ch1H      =  100.0e+3;
const double Nch1      =  200.0;
// GBM 2 Channel
const double ch2L      =  100.0e+3;
const double ch2H      =  300.0e+3;
const double Nch2      =  150.0;
// GBM 3 Channel:
const double ch3L            = 300.0e+3;  //1GeV
const double ch3H      = 1.0e+6; //2 GeV
const double Nch3      =  100.0;
// LAT RANGE
double ch4L      = 1.0e+9; //!MeV this will be change by enph
const double ch4H      = 1.0e+9;  //!GeV 
const double Nch4      =  10.0;

const double ch5L      = 10.0e+9; 
const double ch5H      = cst::enmax;
const double Nch5      =  1.0;
////////////////////////////////////////////////////////////

const double TIME=10.0;
const double EVENTS=100000;
double m_flux[ENERGYSIZE][TIMESIZE];
double timex[TIMESIZE];
double energyx[ENERGYSIZE];

void Burst(int argc, char** argv);

int main(int argc, char** argv)
{
  std::cout<<" ****** GRB and ROOT test ****** "<<std::endl;
  std::cout<<" type test_GRBROOT.exe -help for help on the available options "
	   <<std::endl;
  Burst(argc,argv);
  return 0;
}

void Burst(int argc, char** argv)
{
  // definition and default options 
  bool savef      = false;
  bool save_gif   = false;
  bool save_root  = false;
  bool video_out  = true;
  bool movie      = true;
  bool set_tmax   = false;
  static const char * default_arg="GRBSpectrum";
  int current_arg = 1;
  int Ncycle      = 1;
  int seed        = 0;
  double tmax;
  
  std::vector<double> spectrum,energy,deltae;
  std::string arg_name(default_arg);
  
  while(current_arg < argc)
    {
      arg_name = argv[current_arg]; 
      if ("-help"==arg_name || "-h"==arg_name)
	{
  	  std::cout<<" ** HELP for test_GRBROOT.exe **"<<std::endl;
	  std::cout<<" Description:"<<std::endl;
	  std::cout<<"    GRBROOTTest is a test program for GRB simulation 
studies."<<std::endl;
	  std::cout<<"    This executable uses ROOT to display several 
histograms at run time."<<std::endl;
	  std::cout<<"    It is just an interface to GRBSim, which actually 
does all the job."<<std::endl;
	  std::cout<<" "<<std::endl;
	  std::cout<<" Usage:                                      "<<std::endl;
	  std::cout<<"    test_GRBROOT.exe [options] [option_args] "<<std::endl;
	  std::cout<<" "<<std::endl;
	  std::cout<<" Options:                                    "<<std::endl;
	  std::cout<<"    -help or -h print this help         "<<std::endl;
	  std::cout<<"    -time <time in sec> fix the maximum time "<<std::endl;
	  std::cout<<"    -save Save the parameters for the model and 
\t\t the output in a file (\'GRBdata.txt\')."<<std::endl;
	  std::cout<<"    -root Save the histograms in a in a root file 
\t\t (\'histos.root\')"<<std::endl;
	  std::cout<<"    -seed <seed> set the seed for the random engine"
		   <<std::endl;
	  std::cout<<"    -mute It turns off the graphical output (fastest)"
		   <<std::endl;
	  std::cout<<"    -nomovie It turn off the visualization of the 
\t\t evolution of the flux (faster)"<<std::endl;
	  std::cout<<"    -cycle <N> it runs the program N times 
(with mute option) "<<std::endl;
	  std::cout<<"    -gif  It save the moovie of the evolving flux 
\t\t in a gif file."<<std::endl;
	  std::cout<<" Author:         Nicola Omodei "<<std::endl;
	  exit(0);
	}
      else if("-time"==arg_name)
	{
	  set_tmax=true;
	  tmax = atof(argv[++current_arg]);
	}
      else if("-save"==arg_name)
	{
	  savef=true;
	}
      else if("-root"==arg_name)
	{
	  save_root=true;
	}
      else if("-mute"==arg_name)
	{
	  video_out=false;
	  movie=false;
	}
      else if("-nomovie"==arg_name)
	{
	  movie=false;
	}

      else if("-cycle"==arg_name)
	{
	  video_out=false;
	  Ncycle = atoi(argv[++current_arg]);
	}
      else if("-gif"==arg_name)
	{
	  save_gif=true;
	}
      else if("-seed"==arg_name)
	{
	  seed=atoi(argv[++current_arg]);
	}
      current_arg++;
    }
  
  for(int i=0;i<Ncycle;i++)
    {
      // 1) Initialize ROOT
      std::cout<<"******     Initializing ROOT    ******"<<std::endl;
      static const   TROOT root("GRB", "GRB Simulator");
      // 2) Initialize the GRB simulation, setting the seed
      GRBSim* myGRB = new GRBSim(seed);
      // 3) Build a new GRB
      myGRB->MakeGRB();
      
      // 4) Prepare the time and energy binning 
      if (!set_tmax) tmax=myGRB->Tmax();
      double m_dt=tmax/cst::nstep;
      energy=myGRB->Energy();
      deltae=myGRB->DeltaE();
      ch4L  =myGRB->myParam->EnergyPh();
      TApplication theApp("App",0,0);
      // Logaritmic Binning for energy
      double *e_bins;
      e_bins=(double*)malloc(sizeof(double)*(cst::enstep+1));
  
      for(int i = 0;i<=cst::enstep;i++)
	{
	  e_bins[i]=energy[i]*1.0e-6;
	}
     
      TH1D *h1 = new TH1D("h1","h1",cst::enstep,e_bins);
      
      TCanvas *cc2;
      if(video_out && movie)
	{
	  cc2 = new TCanvas();
	}
      
      double fluence1=0.0;
      double fluence2=0.0;
      double fluence3=0.0;
      double fluence4=0.0;
      double fluence5=0.0;
      double fluenceTOT=0.0;
      std::vector<double> m_time;
      
      // to not have confusion between widows and linux...
      int en,t;
      // 5) STARTS the main loop. For each time, it computs the flux and update the histograms
      for (t=0;t<cst::nstep;t++)
	{
	  m_time.push_back(t*m_dt);
	  spectrum.clear();	  
	  spectrum=myGRB->ComputeFlux(m_time[t]); // is in ph/s/MeV/m^2
	  for (en=0;en<=cst::enstep;en++)
	    {
	      m_flux[en][t]=spectrum[en]*energy[en]/(cst::erg2MeV*1.0e+6);  
	      // m_flux is in erg/s/MeV/m^2
	    }
	  
	  if (video_out && movie) 
	    {
	      for (en=0;en<cst::enstep;en++)
		{
		  h1->SetBinContent(en+1,m_flux[en][t]*1.0e-4); //erg/cm^2/s/MeV
		}
	      cc2->cd();
	      cc2->SetLogx();
	      cc2->SetLogy();
	      
	      h1->SetStats(kFALSE);
	      h1->Draw("AL");
	      h1->GetXaxis()->SetTitle("Energy[MeV]");
	      h1->GetXaxis()->SetTitleSize(0.035);
	      h1->GetXaxis()->SetLabelSize(0.035);
	      h1->GetXaxis()->SetTitleOffset(1.4);
	      h1->GetYaxis()->SetTitle("#nu F#nu[erg/cm^2/s/MeV]");
	      h1->GetYaxis()->SetTitleSize(0.035);
	      h1->GetYaxis()->SetLabelSize(0.035);
	      h1->GetYaxis()->SetTitleOffset(1.4);   
	      h1->SetMaximum(1.0e-2);
	      h1->SetMinimum(1.0e-10);
	      h1->SetLineWidth(4);
	      h1->SetLineColor(4);
	      
	      TLegend *l1 = new TLegend(0.13,0.85,0.85,0.98);
	      l1->SetTextFont(72);
	      l1->SetTextSize(0.05);
	      
	      char labeltime[100];
	      sprintf(labeltime," Time/Tmax = %4g / %4g  ",m_time[t],tmax);
	      l1->AddEntry(h1,labeltime,"");
	      l1->Draw();
	      
	      cc2->Update();
	      char gif[24];
	      if(save_gif) 
		{
		  sprintf(gif, "flux_%.3d.gif", t);
		  cc2->SaveAs(gif);
		}
	    }
	  //myGRB->IFlux returns eV/s/m^2 ; fluences are in erg/cm^2: 
	  fluenceTOT+=myGRB->IFlux(spectrum,cst::enmin,cst::enmax)/(cst::erg2MeV*1.0e+6)*(1.0e-4)*m_dt;
	  fluence1+=myGRB->IFlux(spectrum,ch1L,ch1H)              /(cst::erg2MeV*1.0e+6)*(1.0e-4)*m_dt;
	  fluence2+=myGRB->IFlux(spectrum,ch2L,ch2H)              /(cst::erg2MeV*1.0e+6)*(1.0e-4)*m_dt;
	  fluence3+=myGRB->IFlux(spectrum,ch3L,ch3H)              /(cst::erg2MeV*1.0e+6)*(1.0e-4)*m_dt;
	  fluence4+=myGRB->IFlux(spectrum,ch4L,ch4H)              /(cst::erg2MeV*1.0e+6)*(1.0e-4)*m_dt;
	  fluence5+=myGRB->IFlux(spectrum,ch5L,ch5H)              /(cst::erg2MeV*1.0e+6)*(1.0e-4)*m_dt;
	  // erg/cm^2
	  if (video_out)
	    std::cout<<"Time / tmax = "<<m_time[t]<<"/" <<tmax<<" Ftot [erg]= "
		<<fluenceTOT*(myGRB->Area()*1.0e+4)<<std::endl;
	}
      //////////////////////////////////////////////////
      if(save_gif) 
	{
	  system("gifsicle --delay=5 --loop=5 flux_*.gif > anim.gif");
	  system("rm -f flux_*.gif");
	}
      
      // gSystem->Exec("gifsicle --delay=10 --loop=5 flux_*.gif > anim.gif");
      // gSystem->Exec("rm -f flux_*.gif");
      
      //You can view the animated file anim.gif with Netscape/IE
      //or with gifview as shown below (finish by typing "q" in the window)
      //for more info about gifsicle, gifview, see above url or do
      //gifsicle --help     gifview --help
      //gSystem->Exec("gifview -a anim.gif");
      //////////////////////////////////////////////////

      std::cout<<"**************************************************"<<std::endl;
      if (fluence2!=0) 
	std::cout<<" HardNess Ratio =  "<<fluence3/fluence2<<std::endl;
      std::cout<<"Fluence [erg/cm^2] = "<<fluenceTOT<<std::endl;
      
      std::cout<<"Fluence ("<<ch1L*1.0e-6<<" MeV - "<<ch1H*1.0e-6<<" MeV) [erg/cm^2] = "
	  <<fluence1<<", photons expected = "<<
	fluence1*(cst::erg2MeV)*1.0e+10/(ch1H-ch1L)*6.<<std::endl;
      std::cout<<"Fluence ("<<ch2L*1.0e-6<<" MeV - "<<ch2H*1.0e-6<<" MeV) [erg/cm^2] = "
	  <<fluence2<<", photons expected =  "<<
	fluence2*(cst::erg2MeV)*1.0e+10/(ch2H-ch2L)*6.<<std::endl;
      std::cout<<"Fluence ("<<ch3L*1.0e-6<<" MeV - "<<ch3H*1.0e-6<<" MeV) [erg/cm^2] = "
	  <<fluence3<<", photons expected =  "<<
	fluence3*(cst::erg2MeV)*1.0e+10/(ch3H-ch3L)*6.<<std::endl;
      std::cout<<"Fluence ("<<ch4L*1.0e-6<<" MeV - "<<ch4H*1.0e-6<<" MeV) [erg/cm^2] = "
	  <<fluence4<<", photons expected =  "<<
	fluence4*(cst::erg2MeV)*1.0e+10/(ch4H-ch4L)*6.<<std::endl;
      std::cout<<"Fluence ("<<ch5L*1.0e-6<<" MeV - "<<ch5H*1.0e-6<<" MeV) [erg/cm^2] = "
	  <<fluence5<<", photons expected =  "<<
	fluence5*(cst::erg2MeV)*1.0e+10/(ch5H-ch5L)*6.<<std::endl;

      std::cout<<"Fluence ("<<ch1L*1.0e-6<<" MeV - "<<ch1H*1.0e-6<<" MeV) [ph/cm^2/s] = "
	  <<fluence1*(cst::erg2MeV)*1.0e+6/(ch4H-ch4L)/tmax<<std::endl;
      std::cout<<"Fluence ("<<ch2L*1.0e-6<<" MeV - "<<ch2H*1.0e-6<<" MeV) [ph/cm^2/s] = "
	  <<fluence2*(cst::erg2MeV)*1.0e+6/(ch4H-ch4L)/tmax<<std::endl;
      std::cout<<"Fluence ("<<ch3L*1.0e-6<<" MeV - "<<ch3H*1.0e-6<<" MeV) [ph/cm^2/s] = "
	  <<fluence3*(cst::erg2MeV)*1.0e+6/(ch4H-ch4L)/tmax<<std::endl;
      std::cout<<"Fluence ("<<ch4L*1.0e-6<<" MeV - "<<ch4H*1.0e-6<<" MeV) [ph/cm^2/s] = "
	  <<fluence4*(cst::erg2MeV)*1.0e+6/(ch4H-ch4L)/tmax<<std::endl;
      std::cout<<"Fluence ("<<ch5L*1.0e-6<<" MeV - "<<ch5H*1.0e-6<<" MeV) [ph/cm^2/s] = "
	  <<fluence5*(cst::erg2MeV)*1.0e+6/(ch4H-ch4L)/tmax<<std::endl;
      //////////////////////////////////////////////////
      // Save the results in the paramFile (defined in GRBConstant.h)
      //////////////////////////////////////////////////
      if(savef)
	{
	  myGRB->myParam->Save(!(cst::savef));
	  std::string paramFile = myGRB->myParam->GetParamFile();
	  facilities::Util::expandEnvVar(&paramFile);
	
	  std::ofstream f2(paramFile.c_str(),ios::app);
	  
	  f2<<seed<<std::endl;
	  f2<<tmax<<std::endl;
	  f2<<fluenceTOT*(myGRB->Area()*1.e+4)<<std::endl;
	  f2<<fluenceTOT<<std::endl;
	  f2<<fluence1<<std::endl;
	  f2<<fluence2<<std::endl;
	  f2<<fluence3<<std::endl;
	  f2<<fluence4<<std::endl;
	  f2.close();
	}
      
      //////////////////////////////////////////////////
      //               Plot - Graphics - Histograms   //
      //////////////////////////////////////////////////
      
      if(video_out)
	{
	  for (t=0;t<cst::nstep;t++)
	    {
	      timex[t]=m_time[t];
	    }
	  for (en=0;en<=cst::enstep;en++)
	    {
	      energyx[en]=energy[en];
	    }
	  
	  TCanvas *c1 = new TCanvas("c1","GRB Flux",10,10,700,800);
	  TCanvas *c2 = new TCanvas("c2","GRB Light Curve");
	  TCanvas *c3 = new TCanvas("c3","3D Flux Rapresentation");
	 	  
	  TH1D *fl = new TH1D("fl","fl",cst::enstep,e_bins);
	  TH1D *ph = new TH1D("ph","ph",cst::enstep,e_bins);
	  //TGraph *gfl = new TGraph(cst::enstep);
	  //TGraph *gph = new TGraph(cst::enstep);
	  /*	  
		  TGraph *glct = new TGraph(cst::nstep);
		  TGraph *gch1 = new TGraph(cst::nstep);
		  TGraph *gch2 = new TGraph(cst::nstep);
		  TGraph *gch3 = new TGraph(cst::nstep);
		  TGraph *gch4 = new TGraph(cst::nstep);
	  */
	  TH1D *glct = new TH1D("glct","Ligh Curves",cst::nstep,timex[0],timex[cst::nstep-1]);
	  TH1D *gch1 = new TH1D("gch1","gch1",cst::nstep,timex[0],timex[cst::nstep-1]);
	  TH1D *gch2 = new TH1D("gch2","gch2",cst::nstep,timex[0],timex[cst::nstep-1]);
	  TH1D *gch3 = new TH1D("gch3","gch3",cst::nstep,timex[0],timex[cst::nstep-1]);
	  TH1D *gch4 = new TH1D("gch4","gch4",cst::nstep,timex[0],timex[cst::nstep-1]);
	  TH1D *gch5 = new TH1D("gch5","gch5",cst::nstep,timex[0],timex[cst::nstep-1]);
	  
	  TH2D *h2 = new TH2D("h2","3D flux",cst::enstep,e_bins,cst::nstep,timex[0],timex[cst::nstep-1]);
	  //TH2D *h2 = new TH2D("h2","H2",cst::enstep-1,2,20,cst::nstep-1,0,1);
	  //h1->SetOption("B");
	  
	  double Fv[cst::enstep];
	  double vFv[cst::enstep];
	  double ph_m2[cst::enstep];
	  double lct[cst::nstep];
	  double ch1[cst::nstep];
	  double ch2[cst::nstep];
	  double ch3[cst::nstep];
	  double ch4[cst::nstep];
	  double ch5[cst::nstep];
	  
	  for (t=0;t<cst::nstep;t++)
	    {
	      ch1[t]=0.0;
	      ch2[t]=0.0;
	      ch3[t]=0.0;
	      ch4[t]=0.0;
	      ch5[t]=0.0;
	      lct[t]=0.0;
	    }
	  double loge[cst::enstep];
	  double de;
	  for (en=0;en<cst::enstep;en++)
	    {
	      Fv[en]=0.0;
	      vFv[en]=0.0;
	      ph_m2[en]=0.0;
	      loge[en]=log10(energyx[en]);
	      de=(energyx[en+1]-energyx[en]);

	      for (t=0;t<cst::nstep;t++) 
		{
		  // m_flux is in erg/s/MeV/m^2
		  Fv[en]+=(m_flux[en][t])*m_dt;  // is in erg/m^2/MeV
		  vFv[en]+=(m_flux[en][t])*(energyx[en])*1.0e-4;
		  // is in erg/s/cm^2
		  // ph_m2[en]+=(m_flux[en][t])*(cst::erg2MeV*1.0e+6)/
		  //energyx[en]; 
		  // is in ph/MeV/m^2/s
		  ph_m2[en]+=(m_flux[en][t])*(cst::erg2MeV*1.0e-4)
		    /(energyx[en]*1.0e-3);
		  // is in ph/KeV/cm^2/s
		  lct[t]+=(m_flux[en][t])*(cst::erg2MeV)
		    /energyx[en]*de*m_dt;
		  //ph/m^2
		  
		  if (energyx[en]<ch1H && energyx[en]>=ch1L)
		    {
		      ch1[t]+=(m_flux[en][t])*(cst::erg2MeV)
			/energyx[en]*de*m_dt;
		      // is in ph/m^2
		    }
		  if (energyx[en]<ch2H &&energyx[en]>=ch2L)
		    {
		      ch2[t]+=(m_flux[en][t])*(cst::erg2MeV)
			/energyx[en]*de*m_dt;
		      // is in ph/m^2
		    }
		  if (energyx[en]<ch3H &&energyx[en]>=ch3L)
		    {
		      ch3[t]+=(m_flux[en][t])*(cst::erg2MeV)
			/energyx[en]*de*m_dt;
		      // is in ph/m^2

 		    }
		  if (energyx[en]<ch4H &&energyx[en]>=ch4L)
		    {
		      ch4[t]+=(m_flux[en][t])*(cst::erg2MeV) 
			/energyx[en]*de*m_dt;
		      // is in ph/m^2
		    }
		  if (energyx[en]<ch5H &&energyx[en]>=ch5L)
		    {
		      ch5[t]+=(m_flux[en][t])*(cst::erg2MeV) 
			/energyx[en]*de*m_dt;
		      // is in ph/m^2
		    }
		  h2->SetBinContent(en+1,t+1,m_flux[en][t]);
		}
	    }
	  for (en=0;en<cst::enstep;en++)
	    {
	      //	      fl->SetBinContent(en+1,Fv[en]); 
	      fl->SetBinContent(en+1,vFv[en]); 
	      ph->SetBinContent(en+1,ph_m2[en]);
	    }
	  TRandom *noise = new TRandom();
	  for (t=0;t<cst::nstep;t++)
	    {
	      glct->SetBinContent(t+1,lct[t]);
	      gch1->SetBinContent(t+1,ch1[t]+noise->Poisson(Nch1));
	      gch2->SetBinContent(t+1,ch2[t]+noise->Poisson(Nch2));
	      gch3->SetBinContent(t+1,ch3[t]+noise->Poisson(Nch3));
	      gch4->SetBinContent(t+1,ch4[t]+noise->Poisson(Nch4));
	      gch5->SetBinContent(t+1,ch5[t]+noise->Poisson(Nch5));
	    }
	  
	  gch1->SetLineColor(2);
	  gch2->SetLineColor(3);
	  gch3->SetLineColor(4);
	  gch4->SetLineColor(8);
	  gch5->SetLineColor(6);
	  /*
	    gch1->SetFillColor(2);
	    gch2->SetFillColor(3);
	    gch3->SetFillColor(4);
	    gch4->SetFillColor(6);
	    
	    gch1->SetFillStyle(4050);
	    gch2->SetFillStyle(4050);
	    gch3->SetFillStyle(4050);
	    gch4->SetFillStyle(4050);
	  */
	  c1->Divide(1,2);

	  TPad *c1_1 = (TPad*) c1->GetPrimitive("c1_1");
	  TPad *c1_2 = (TPad*) c1->GetPrimitive("c1_2");
	
	  c1->cd(1);
	  c1_1->SetLogx();
	  c1_1->SetLogy();  
	  c1_1->SetGridx();
	  c1_1->SetGridy();
	  ph->SetLineColor(4);
	  ph->SetLineWidth(4);
	  ph->SetStats(kFALSE);
	  ph->GetXaxis()->SetTitle("Energy[MeV]");
	  ph->GetXaxis()->SetTitleSize(0.035);
	  ph->GetXaxis()->SetLabelSize(0.035);
	  ph->GetXaxis()->SetTitleOffset(1.4);
	  ph->GetYaxis()->SetTitle("Photons/cm^2/KeV/s]");
	  ph->GetYaxis()->SetTitleSize(0.035);
	  ph->GetYaxis()->SetLabelSize(0.05);
	  ph->GetYaxis()->SetTitleOffset(1.4);   
	  ph->SetMinimum(1.0e-10);
	  ph->Draw();
	  
	  TF1 *fit1= new TF1("fit1","[0]*x^[1]",cst::enmin*1.0e-6,cst::enmax*1.0e-6);
	  fit1->SetParameters(1,-2);
	  fit1->FixParameter(1,-2);
	  fit1->SetLineColor(3);
	  ph->Fit("fit1","R+");
	  ph->Draw("AL");
	  
	  c1->cd(2);
	  glct->SetLineColor(5);
	  glct->SetLineStyle(3);
	  
	  c1_2->SetLogx();
	  c1_2->SetLogy();
	  c1_2->SetGridx();
	  c1_2->SetGridy();
	  //c1_1->SetLogx();
	  //c1_1->SetLogy();
	  fl->SetLineColor(4);
	  fl->SetLineWidth(4);
	  fl->SetStats(kFALSE);
	  fl->GetXaxis()->SetTitle("Energy[MeV]");
	  fl->GetXaxis()->SetTitleSize(0.035);
	  fl->GetXaxis()->SetLabelSize(0.035);
	  fl->GetXaxis()->SetTitleOffset(1.4);
	  fl->GetYaxis()->SetTitle("#nu F#nu[erg/s/cm^2]");  
	  	  
	  fl->GetYaxis()->SetTitleSize(0.035);
	  fl->GetYaxis()->SetLabelSize(0.05);
	  fl->GetYaxis()->SetTitleOffset(1.4);   	  
	  fl->Draw("AL");
	  
	  c1->Update();
	  //c1b->Update();
	  char legendName0[100];
	  char legendName1[100];
	  char legendName2[100];
	  char legendName3[100];
	  char legendName4[100];
	  char legendName5[100];
	  //char legendName5[100];
	  
	  c2->cd(); 
	  c2->SetFillColor(1);
	  c2->SetFrameLineColor(5);
	  c2->SetFrameLineWidth(2);
	  c2->SetBorderSize(2);
	  c2->SetLeftMargin(0.12);
	  c2->SetTopMargin(0.12);
	  c2->SetBottomMargin(0.12);
	  glct->SetStats(kFALSE);
	  glct->GetXaxis()->SetAxisColor(5);
	  glct->GetXaxis()->SetLabelColor(5);
	  glct->GetXaxis()->SetTitleColor(5);
	  glct->GetYaxis()->SetAxisColor(5);
	  glct->GetYaxis()->SetLabelColor(5);
	  glct->GetYaxis()->SetTitleColor(5);
	  glct->GetYaxis()->SetAxisColor(5);
	  
	  glct->GetXaxis()->SetTitle("Time [sec]");
	  glct->GetXaxis()->SetLabelSize(0.035);
	  glct->GetXaxis()->SetTitleSize(0.035);
	  glct->GetXaxis()->SetTitleOffset(1.4);
	  glct->GetYaxis()->SetTitle("Counts [ph/m^2]");
	  glct->GetYaxis()->SetLabelSize(0.035);
	  glct->GetYaxis()->SetTitleOffset(1.2);
	  glct->Draw("AL");
	  	  
	  gch1->Draw("same");
	  gch2->Draw("same");
	  gch3->Draw("same");
	  gch4->Draw("same");
	  gch5->Draw("same");
	  
	  //   LEGEND
	  TLegend *leg = new TLegend(0.68,0.80,0.98,0.98);
	  leg->SetTextFont(11);
	  sprintf(legendName0,"(%6g GeV - %6g GeV) ",cst::enmin*1.0e-9,cst::enmax*1.0e-9);
	  sprintf(legendName1,"(%6g GeV - %6g GeV) ",ch1L*1.0e-9,ch1H*1.0e-9);
	  sprintf(legendName2,"(%6g GeV - %6g GeV) ",ch2L*1.0e-9,ch2H*1.0e-9);
	  sprintf(legendName3,"(%6g GeV - %6g GeV) ",ch3L*1.0e-9,ch3H*1.0e-9);
	  sprintf(legendName4,"(%6g GeV - %6g GeV) ",ch4L*1.0e-9,ch4H*1.0e-9);
	  sprintf(legendName5,"(%6g GeV - %6g GeV) ",ch5L*1.0e-9,ch5H*1.0e-9);
	  //  sprintf(legendName5,"(%6g GeV - %6g GeV) ",ch1L,ch1H);
	  leg->AddEntry(glct,legendName0,"l");
	  leg->AddEntry(gch1,legendName1,"l");
	  leg->AddEntry(gch2,legendName2,"l");
	  leg->AddEntry(gch3,legendName3,"l");
	  leg->AddEntry(gch4,legendName4,"l");
	  leg->AddEntry(gch5,legendName5,"l");
	  leg->Draw();
	  
	  
	  c2->Update();
	  
	  c3->cd();
	  c3->SetLogx();
	  c3->SetLogz(); 
	  h2->GetXaxis()->SetTitle("Energy [MeV]");
	  h2->GetXaxis()->SetTitleOffset(2);
	  h2->GetYaxis()->SetTitle("Time[sec]");
	  h2->GetYaxis()->SetLabelOffset(0);
	  h2->GetYaxis()->SetTitleOffset(2);
	  h2->GetZaxis()->SetTitle("vFv[erg/s/m^2/MeV]");
	  h2->GetZaxis()->SetLabelSize(0.05);
	  h2->GetZaxis()->SetTitleOffset(1.3);
	  h2->GetXaxis()->CenterTitle(1);
	  h2->GetYaxis()->CenterTitle(1);
	  h2->GetZaxis()->CenterTitle(1);
	  h2->Draw("surf3");
	  c3->Update();
	  /*
	    c4->cd();  
	    c4->SetLogx();
	    c4->SetLogz();
	    h2->Draw("CONT");
	    c4->Update();
	  */
	  if(save_root)
	    {
	      std::cout<<"SAVING The Histograms in a root files..."<<std::endl;
	      TFile f1("histos.root","recreate");
	      //LigtCourves
	      glct->Write();
	      gch1->Write();
	      gch2->Write();
	      gch3->Write();
	      gch4->Write();
	      gch5->Write();
	      //3d Plot
	      h2->Write(); 
	      //Flux
	      ph->Write(); 
	      fl->Write();
	      c1->Write();
	      c2->Write();
	      c3->Write();
	      //
	      f1.Close();
	      
	    }
	}
      ////////////////////////////////////////////
      
      if (video_out) theApp.Run();
      delete myGRB;
      //To kill the application press control C //
    }
}
