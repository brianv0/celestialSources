/*!\class GRBROOTTest
 * \brief Test program for GRB simulation studies.
 * 
 * This executable uses ROOT to display several histograms at run time.
 * To use it, type on the prompt ./test_GRBROOT.exe
 * After a period of initialization, a canvas pops up showing the 
 * complete evolution of the burst in function of the time. 
 * At the end, 4 new canvas pops up, showing several summary 
 * histograms.
 */
#include <iterator>
#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>
#include "facilities/Util.h"
//Include files for spectrum...

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
//#include "BranchObject.h"
#include "TFile.h" 

//GRB include files...
#include "../../GRB/GRBShell.h"
#include "../../GRB/GRBShock.h"
#include "../../GRB/GRBConstants.h"
#include "../../GRB/GRBSim.h"
//#include "../../src/GRB/GRBSpectrum.h"
#include "CLHEP/Random/RandFlat.h"
using namespace std;

//! The size of the energy array is taken from GRBConstant.h file in FluxSvc
#define ENERGYSIZE cst::enstep
//! The size of the time array is taken from GRBConstant.h file in FluxSvc
#define TIMESIZE cst::nstep
////////////////////////////////////////////////////////////
// Channel for the output Light Curves
// BATSE 1 Channel:
const double ch1L      = 20.0e+3;
const double ch1H      = 50.0e+3;
// BATSE 2 Channel
const double ch2L      = 50.0e+3;
const double ch2H      = 100.0e+3;
// BATSE 3 Channel:
const double ch3L      = 100.0e+3;
const double ch3H      = 300.0e+3;
// LAT RANGE
double ch4L      = 10.0e+6; //!MeV
const double ch4H      = cst::enmax;
// GLAST Energies ?
//const double ch1L      = 0.001e+9;
//const double ch1H      = 0.0025e+9;
//const double ch2L      = 0.0025e+9;
//const double ch2H      = 0.005e+9;
//const double ch3L      = 0.005e+9;
//const double ch3H      = 0.01e+9;
//const double ch4L      = 0.01e+9;
//const double ch4H      = 1.00e+9;

// GBM Energy Range:
//const double ch3L      = 25.0e+3;
//const double ch3H      = 10.0e+6;
// LAT Energy Range
//const double ch4L      = 20.0e+6;
//const double ch4H      = 300.0e+9;
////////////////////////////////////////////////////////////

const double TIME=10.0;
const double EVENTS=100000;
double m_flux[ENERGYSIZE][TIMESIZE];
double timex[TIMESIZE];
double energyx[ENERGYSIZE];

//! It starts the simulation of a new Gamma-Ray Burst
void Burst(int argc, char** argv);
int main(int argc, char** argv)
{
  Burst(argc,argv);
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

void Burst(int argc, char** argv)
{
  // definition and default options // 
  bool savef      = false;
  bool save_gif   = false;
  bool video_out  = true;
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
      if("-time"==arg_name)
	{
	  set_tmax=true;
	  tmax = atof(argv[++current_arg]);
	}
      else if("-save"==arg_name)
	{
	  savef=true;
	}
      else if("-mute"==arg_name)
	{
	  video_out=false;
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
      GRBSim* _myGRB = new GRBSim(seed);
      /// This initializes the GRB simulation...
      cout<<"******     Initializing ROOT    ******"<<endl;
      static const   TROOT root("GRB", "GRB Simulator");
      //TROOT->Stat()
      _myGRB->MakeGRB();
      
      if (!set_tmax) tmax=_myGRB->Tmax();
      
      cout<<" MAX TIME = "<<tmax<<endl;
      double m_dt=tmax/cst::nstep;
      
      energy=_myGRB->Energy();
      deltae=_myGRB->DeltaE();
      ch4L  =_myGRB->myParam->EnergyPh();
      TApplication theApp("App",0,0);
      
      
      //------------------------------------------------------//
      ////////////////  Compute the Total Flux /////////////////
      //------------------------------------------------------//
      // Logaritmic Binning
      //      double e_delta=pow(cst::enmax/cst::enmin,1.0/cst::enstep);
      double *e_bins;
      e_bins=(double*)malloc(sizeof(double)*(cst::enstep+1));
  
      for(int i = 0;i<=cst::enstep;i++)
	{
	  //	  e_bins[i]=cst::enmin*pow(e_delta,1.0*i)*1.0e-6;
	  e_bins[i]=energy[i]*1.0e-6;
	}
      
      TH1D *h1 = new TH1D("h1","h1",cst::enstep,e_bins);//log10(cst::enmin)-6,log10(cst::enmax)-6);
      
      TCanvas *cc2;
      if(video_out)
	{
	  cc2 = new TCanvas();
	 	  //	  h1->GetXaxis()->SetTitle("Energy [MeV]");
	  //  h1->GetXaxis()->SetLabelSize(0.05);
	  // h1->GetXaxis()->SetTitleOffset(1.5);
	  //  h1->GetYaxis()->SetTitle("flux[ph/s/eV/m^2]");  
	  //  h1->GetYaxis()->SetTitle("F#nu[erg/s/cm^2]");  
	  //  h1->SetFillColor(5);
	}
      
      double fluence1=0.0;
      double fluence2=0.0;
      double fluence3=0.0;
      double fluence4=0.0;
      double fluenceTOT=0.0;
      std::vector<double> m_time;
      
      for (int t=0;t<cst::nstep;t++)
	{
	  m_time.push_back(t*m_dt);
	  // Compute the flux @ time
	  spectrum.clear();	  
	  spectrum=_myGRB->ComputeFlux(m_time[t]); // is in ph/s/MeV/m^2
	  for (int en=0;en<=cst::enstep;en++)
	    {
	      m_flux[en][t]=spectrum[en]*energy[en]/(cst::erg2MeV*1.0e+6);  
	      // m_flux is in erg/s/MeV/m^2
	    }
	  
	  if (video_out) 
	    {
	      for (int en=0;en<cst::enstep;en++)
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
	  //_myGRB->IFlux returns eV/s/m^2 ; fluences are in erg/cm^2: 
	  fluenceTOT+=_myGRB->IFlux(spectrum,cst::enmin,cst::enmax)/(cst::erg2MeV*1.0e+6)*(1.0e-4)*m_dt;
	  fluence1+=_myGRB->IFlux(spectrum,ch1L,ch1H)              /(cst::erg2MeV*1.0e+6)*(1.0e-4)*m_dt;
	  fluence2+=_myGRB->IFlux(spectrum,ch2L,ch2H)              /(cst::erg2MeV*1.0e+6)*(1.0e-4)*m_dt;
	  fluence3+=_myGRB->IFlux(spectrum,ch3L,ch3H)              /(cst::erg2MeV*1.0e+6)*(1.0e-4)*m_dt;
	  fluence4+=_myGRB->IFlux(spectrum,ch4L,ch4H)              /(cst::erg2MeV*1.0e+6)*(1.0e-4)*m_dt;
	  // erg/cm^2
	  if (video_out) 
	    cout<<"Time / tmax = "<<m_time[t]<<"/" <<tmax<<" Ftot [erg]= "
		<<fluenceTOT*(_myGRB->Area()*1.0e+4)<<endl;
	}
      //////////////////////////////////////////////////
      if(save_gif) 
	{
	  system("gifsicle --delay=5 --loop=5 flux_*.gif > anim.gif");
	  system("rm -f flux_*.gif");
	}
      
      //      gSystem->Exec("gifsicle --delay=10 --loop=5 flux_*.gif > anim.gif");
      // gSystem->Exec("rm -f flux_*.gif");
      
      //You can view the animated file anim.gif with Netscape/IE
      //or with gifview as shown below (finish by typing "q" in the window)
      //for more info about gifsicle, gifview, see above url or do
      // gifsicle --help     gifview --help
      //  gSystem->Exec("gifview -a anim.gif");
      //////////////////////////////////////////////////

      cout<<"**************************************************"<<endl;
      if (fluence2!=0) 
	cout<<" HardNess Ratio =  "<<fluence3/fluence2<<endl;
      cout<<"Fluence [erg/cm^2] = "<<fluenceTOT<<endl;
      
      cout<<"Fluence ("<<ch1L*1.0e-6<<" MeV - "<<ch1H*1.0e-6<<" MeV) [erg/cm^2] = "
	  <<fluence1<<", photons expected = "<<
	fluence1*(cst::erg2MeV)*1.0e+10/(ch1H-ch1L)*6.<<endl;
      cout<<"Fluence ("<<ch2L*1.0e-6<<" MeV - "<<ch2H*1.0e-6<<" MeV) [erg/cm^2] = "
	  <<fluence2<<", photons expected =  "<<
	fluence2*(cst::erg2MeV)*1.0e+10/(ch2H-ch2L)*6.<<endl;
      cout<<"Fluence ("<<ch3L*1.0e-6<<" MeV - "<<ch3H*1.0e-6<<" MeV) [erg/cm^2] = "
	  <<fluence3<<", photons expected =  "<<
	fluence3*(cst::erg2MeV)*1.0e+10/(ch3H-ch3L)*6.<<endl;
      cout<<"Fluence ("<<ch4L*1.0e-6<<" MeV - "<<ch4H*1.0e-6<<" MeV) [erg/cm^2] = "
	  <<fluence4<<", photons expected =  "<<
	fluence4*(cst::erg2MeV)*1.0e+10/(ch4H-ch4L)*6.<<endl;
      //////////////////////////////////////////////////
      // Save the results in the GRBlogfile...        //
      //////////////////////////////////////////////////
      if(savef)
	{
	  _myGRB->myParam->Save(!(cst::savef));
	  std::string paramFile = _myGRB->myParam->GetParamFile();
	  //"$(GRBROOT)/src/test/GRBlog.txt";
	  //	  std::string paramFile2 = "$(GRBROOT)/src/test/GRBdata.txt";
	  
	  //facilities::Util::expandEnvVar(&paramFile);
	  facilities::Util::expandEnvVar(&paramFile);
	  
	  //ofstream f1(paramFile.c_str(),ios::app);
	  ofstream f2(paramFile.c_str(),ios::app);
	  /*
	  if (! f1.is_open()) 
	    {
	      cout<<"Error Opening $(GRBROOT)/src/test/GRBlog.txt\n";
	      //TODO LIST: still need to remove this exit, without getting a core dump!
	      exit(1);
	    }
	  cout<<"Save results into the file: "<<paramFile.c_str()<<endl;
	  cout<<"*******************************************"<<endl;
	  
	  f1<<"Tmax = "<<tmax<<" Ftot [erg]= "<<fluence1<<endl;
	  f1<<"Fluence [erg/cm^2] = "<<fluence1/_myGRB->Area()<<endl;
	  f1<<"Fluence ("<<ch1L*1.0e-6<<" MeV - "<<ch1H*1.0e-6<<" MeV) [erg/cm^2] = "<<fluence2/_myGRB->Area()<<endl;
	  f1<<"Fluence ("<<ch2L*1.0e-6<<" MeV - "<<ch2H*1.0e-6<<" MeV) [erg/cm^2] = "<<fluence3/_myGRB->Area()<<endl;
	  f1<<"Fluence ("<<ch3L*1.0e-6<<" MeV - "<<ch3H*1.0e-6<<" MeV) [erg/cm^2] = "<<fluence4/_myGRB->Area()<<endl;
	  f1<<"Fluence ("<<ch4L*1.0e-6<<" MeV - "<<ch4H*1.0e-6<<" MeV) [erg/cm^2] = "<<fluence5/_myGRB->Area()<<endl;
	  f1.close();
	  */
	  f2<<seed<<endl;
	  f2<<tmax<<endl;
	  f2<<fluenceTOT*(_myGRB->Area()*1.e+4)<<endl;
	  f2<<fluenceTOT<<endl;
	  f2<<fluence1<<endl;
	  f2<<fluence2<<endl;
	  f2<<fluence3<<endl;
	  f2<<fluence4<<endl;
	  f2.close();
	}
      
      //////////////////////////////////////////////////
      //               Plot - Graphics - Histograms   //
      //////////////////////////////////////////////////
      
      if(video_out)
	{
	  for (int t=0;t<cst::nstep;t++)
	    {
	      timex[t]=m_time[t];
	    }
	  for (int en=0;en<=cst::enstep;en++)
	    {
	      energyx[en]=energy[en];
	    }
	  
	  TCanvas *c1 = new TCanvas("c1","GRB Flux",10,10,700,800);
	  //TCanvas *c1b = new TCanvas("c1b","Photons Count");
	  TCanvas *c2 = new TCanvas("c2","GRB Light Curve");
	  //TCanvas *c4 = new TCanvas("c4","Countour Plot");
	  TCanvas *c5 = new TCanvas("c5","3D Flux Rapresentation");
	  //cout<<"xmin = "<< energyx[0] << "   xmax = "<<energyx[cst::nstep-1]<< "  ymin = "<<timex[0] <<"   ymax =" << timex[cst::nstep-1] <<endl;
	 	  
	  TH1D *fl = new TH1D("fl","fl",cst::enstep,e_bins);
	  TH1D *ph = new TH1D("ph","ph",cst::enstep,e_bins);
	  //TGraph *gfl = new TGraph(cst::enstep);
	  //TGraph *gph = new TGraph(cst::enstep);
	  
	  TGraph *glct = new TGraph(cst::nstep);
	  TGraph *gch1 = new TGraph(cst::nstep);
	  TGraph *gch2 = new TGraph(cst::nstep);
	  TGraph *gch3 = new TGraph(cst::nstep);
	  TGraph *gch4 = new TGraph(cst::nstep);
	  
	  TH2D *h2 = new TH2D("h2","H2",cst::enstep,e_bins,cst::nstep,timex[0],timex[cst::nstep-1]);
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
	  
	  for (int t=0;t<cst::nstep;t++)
	    {
	      ch1[t]=0.0;
	      ch2[t]=0.0;
	      ch3[t]=0.0;
	      ch4[t]=0.0;
	      lct[t]=0.0;
	    }
	  double loge[cst::enstep];
	  double de;
	  for (int en=0;en<cst::enstep;en++)
	    {
	      Fv[en]=0.0;
	      vFv[en]=0.0;
	      ph_m2[en]=0.0;
	      loge[en]=log10(energyx[en]);
	      de=(energyx[en+1]-energyx[en]);

	      for (int t=0;t<cst::nstep;t++) 
		{
		  // m_flux is in erg/s/MeV/m^2
		  Fv[en]+=(m_flux[en][t])*m_dt;  // is in erg/m^2/MeV
		  //		  vFv[en]+=(m_flux[en][t])*m_dt*(energyx[en])*1.0e-6; // is in erg/m^2
		  vFv[en]+=(m_flux[en][t])*(energyx[en])*cst::erg2MeV*1.0e-7; // is in KeV/s/cm^2
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
		  h2->SetBinContent(en+1,t+1,m_flux[en][t]);
		}
	    }
	  for (int en=0;en<cst::enstep;en++)
	    {
	      //	      fl->SetBinContent(en+1,Fv[en]); 
	      fl->SetBinContent(en+1,vFv[en]); 
	      ph->SetBinContent(en+1,ph_m2[en]);
	    }
	  //	  gfl=new TGraph(cst::enstep,energyx,vFv);
	  //  gph=new TGraph(cst::enstep,energyx,ph_m2);
	  
	  //glc=new TGraph(cst::nstep,timex,lct);
	  glct = new TGraph(cst::nstep,timex,lct);
	  gch1 = new TGraph(cst::nstep,timex,ch1);
	  gch2 = new TGraph(cst::nstep,timex,ch2);
	  gch3 = new TGraph(cst::nstep,timex,ch3);
	  gch4 = new TGraph(cst::nstep,timex,ch4);
	  
	  //	  gfl->SetLineWidth(2);
	  // gph->SetLineWidth(2);
	  //glc->SetLineWidth(2);
	  // gfl->SetLineColor(4);
	  // gph->SetLineColor(4);
   
	  gch1->SetLineColor(2);
	  gch2->SetLineColor(3);
	  gch3->SetLineColor(4);
	  gch4->SetLineColor(5);
	  
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
	  //	  fl->GetYaxis()->SetTitle("F#nu[erg/m^2/MeV]");  
	  fl->GetYaxis()->SetTitle("#nu F#nu[keV/s/cm^2]");  
	  //	  fl->GetYaxis()->SetTitle("#nu F#nu[erg/m^2]");  
	  
	  fl->GetYaxis()->SetTitleSize(0.035);
	  fl->GetYaxis()->SetLabelSize(0.05);
	  fl->GetYaxis()->SetTitleOffset(1.4);   	  
	  fl->Draw("AL");
	  
	  c1->Update();
	  //c1b->Update();
	  
	  char legendName1[100];
	  char legendName2[100];
	  char legendName3[100];
	  char legendName4[100];
	  //char legendName5[100];
	  
	  c2->cd();
	  glct->Draw("ALP");
	  glct->SetTitle("Light Curves");
	  glct->GetXaxis()->SetTitle("Time [sec]");
	  glct->GetXaxis()->SetTitleSize(0.035);
	  glct->GetXaxis()->SetLabelSize(0.035);
	  glct->GetXaxis()->SetTitleOffset(1.4);
	  glct->GetYaxis()->SetTitle("Counts [ph/m^2]");
	  glct->GetYaxis()->SetTitleSize(0.035);
	  glct->GetYaxis()->SetLabelSize(0.035);
	  glct->GetYaxis()->SetTitleOffset(1.4);
	  glct->Draw("AP");
	  
	  gch1->Draw("LP");
	  gch2->Draw("LP");
	  gch3->Draw("LP");
	  gch4->Draw("LP");
	  
	  
	  TLegend *leg = new TLegend(0.6,0.75,0.98,0.98);
	  leg->SetTextFont(11);
	  
	  sprintf(legendName1,"(%6g GeV - %6g GeV) ",ch1L*1.0e-9,ch1H*1.0e-9);
	  sprintf(legendName2,"(%6g GeV - %6g GeV) ",ch2L*1.0e-9,ch2H*1.0e-9);
	  sprintf(legendName3,"(%6g GeV - %6g GeV) ",ch3L*1.0e-9,ch3H*1.0e-9);
	  sprintf(legendName4,"(%6g GeV - %6g GeV) ",ch4L*1.0e-9,ch4H*1.0e-9);
	  //  sprintf(legendName5,"(%6g GeV - %6g GeV) ",ch1L,ch1H);
	  
	  leg->AddEntry(gch1,legendName1,"l");
	  leg->AddEntry(gch2,legendName2,"l");
	  leg->AddEntry(gch3,legendName3,"l");
	  leg->AddEntry(gch4,legendName4,"l");
	  //  leg->AddEntry(gch5,legendName5,"l");
	  //leg->AddEntry(glc,"Sum","l");
	  leg->Draw();
	  
	  
	  c2->Update();
	  
	  c5->cd();
	  c5->SetLogx();
	  c5->SetLogz();
	  h2->Draw("surf");
	  h2->GetXaxis()->SetTitle("Energy [MeV]");
	  //  h2->GetYaxis()->SetTitleSize(0.05);
	  //  h2->GetXaxis()->SetLabelSize(0.05);
	  h2->GetXaxis()->SetTitleOffset(1.2);
	  h2->GetYaxis()->SetTitle("Time[sec]");
	  //h2->GetYaxis()->SetTitleSize(0.035);
	  //  h2->GetYaxis()->SetLabelSize(0.05);
	  h2->GetYaxis()->SetTitleOffset(1.2);
	  h2->GetZaxis()->SetTitle("vFv[erg/s/m^2/MeV]");
	  // h2->GetYaxis()->SetTitleSize(0.035);
	  h2->GetZaxis()->SetLabelSize(0.05);
	  h2->GetZaxis()->SetTitleOffset(1.3);
	  h2->Draw("surf3");
	  c5->Update();
	  /*
	    c4->cd();  
	    c4->SetLogx();
	    c4->SetLogz();
	    h2->Draw("CONT");
	    c4->Update();
	  */

	}
      ////////////////////////////////////////////
      
      cout<<"YES"<<endl;
      if (video_out) theApp.Run();
      delete _myGRB;
      //To kill the application press control C //
    }
}
