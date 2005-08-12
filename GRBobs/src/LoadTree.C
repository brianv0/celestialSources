
void LoadTree(int type=0, double tmin=0, double tmax=0, TString cuts="")
{
 gStyle->SetCanvasColor(10);
 gStyle->SetPalette(1,0);
 TFile *f = new TFile("GRBTree.root");
 TTree *t = (TTree*) f->Get("Events");
 t->SetDirectory(0);
 //////////////////////////////////////////////////
 TCanvas *LAT_DISPLAY_MC = new TCanvas("LAT_DISPLAY_MC","LAT DISPLAY_MC",800,800);
 // LAT_DISPLAY_MC->EditorBar();
 LAT_DISPLAY_MC->Divide(1,2);
 TPad *LAT_DISPLAY_MC_UP = LAT_DISPLAY_MC->cd(1);
 TPad *LAT_DISPLAY_MC_DO = LAT_DISPLAY_MC->cd(2);
 LAT_DISPLAY_MC_UP->Divide(2,1);
 LAT_DISPLAY_MC_UP->cd(1);
 
 tmin = (tmin==0) ?  t->GetMinimum("Time") : tmin ;
 tmax = (tmax==0) ?  t->GetMaximum("Time") : tmax ;
 
 double emin = 30.0;
 double emax = 300000.0;
 int Ebin    = 20;
 double de = pow(emax/emin,1./Ebin);
 double *e = new double[Ebin +1];
 for(int i = 0; i<=Ebin; i++)
   {
     e[i] = emin*pow(de,1.0*i); //keV
   }
 
 LightCurve = new TH1D("LightCurve","LightCurve",1000,tmin,tmax);
 Spectrum   = new TH1D("Spectrum","Spectrum",Ebin,e);
 
 t->Draw("Time>>LightCurve",cuts);
 LAT_DISPLAY_MC_UP->cd(2);
 gPad->SetLogx();
 gPad->SetLogy();
 t->Draw("Energy>>Spectrum",cuts,"E");
 for(int ei = 1; ei<=Ebin; ei++)
   {
     double ene = Spectrum->GetBinCenter(ei);
     double ne  = Spectrum->GetBinContent(ei);
     double Ene = Spectrum->GetBinError(ei);
     Spectrum->SetBinContent(ei,pow(ene,type) * ne);
     Spectrum->SetBinContent(ei,pow(ene,type) * ne);
     Spectrum->SetBinError(ei,pow(ene,type) * Ene); //ph/cm^2/MeV
 
   }
 Spectrum->SetMarkerStyle(4);
 Spectrum->Draw("e");
 
 LAT_DISPLAY_MC_DO->cd();
 // LAT_DISPLAY_MC_DO->SetLogz();
 // SM = new TH2D("SkyMap","SkyMap",360,-180,180,180,-90,90);
 // t->Draw("B:L>>SkyMap(360,-180,180,180,-90,90)","","AITOFF");
 t->Draw("sin(Theta/2)*cos(Phi):sin(Theta/2)*sin(Phi)>>SkyMap(100,-1,1,100,-1,1)",cuts);
 
 double conv  =  TMath::Pi()/360.0;
 TEllipse *fov10 = new TEllipse(0.,0.,sin( 10.* conv));
 TEllipse *fov20 = new TEllipse(0.,0.,sin( 20. * conv));
 TEllipse *fov30 = new TEllipse(0.,0.,sin( 30. * conv));
 TEllipse *fov45 = new TEllipse(0.,0.,sin( 45. * conv));
 TEllipse *fov90 = new TEllipse(0.,0.,sin( 90. * conv));
 TEllipse *fov180 = new TEllipse(0.,0.,sin( 180. * conv));

 fov10->SetLineStyle(3);
 fov20->SetLineStyle(3);
 fov30->SetLineStyle(3);
 fov45->SetLineStyle(3);
 fov90->SetLineStyle(4);
 fov180->SetLineStyle(1);

 fov10->SetLineColor(2);
 fov20->SetLineColor(2);
 fov30->SetLineColor(2);
 fov45->SetLineColor(2);
 fov90->SetLineColor(2);
 fov180->SetLineColor(2);

 fov10->Draw();
 fov20->Draw();
 fov30->Draw();
 fov45->Draw();
 fov90->Draw();
 fov180->Draw();

 LAT_DISPLAY_MC->cd();
}
