{
  gStyle->SetOptStat(0);
  TCut cut0="Fssc_Fsyn<=0 && theta>30 && beta>-3.0";
  TCut cut1=cut0 +" detected>=1";
  TCut cut2=cut0 +" detected>=5";
  TCut cut3=cut0 +" detected>=10";
  
  gStyle->SetCanvasColor(10);

  TFile *f = new TFile("GRBCatalog.root","OPEN");
  TTree *t = (TTree*) gROOT->FindObject("GrbCatalog");
  new TCanvas("c1","c1");



  t->Draw("beta>>b1",cut0,"");       
  t->Draw("beta>>b2",cut1,"same");       
  t->Draw("beta>>b3",cut2,"same");       
  t->Draw("beta>>b4",cut3,"same");  
  TH1D *b1 = (TH1D*)gROOT->FindObject("b1");
  TH1D *b2 = (TH1D*)gROOT->FindObject("b2");
  TH1D *b3 = (TH1D*)gROOT->FindObject("b3");
  TH1D *b4 = (TH1D*)gROOT->FindObject("b4");
  b1->SetLineColor(1);
  b1->SetLineStyle(2);
  b2->SetLineColor(2);
  b3->SetLineColor(3);
  b4->SetLineColor(4);
  TLegend *leg = new TLegend(0.6,0.7,0.99,0.99,"Applied cuts","brNDC");
  leg->AddEntry(b1,cut0);
  leg->AddEntry(b2,cut1);
  leg->AddEntry(b3,cut2);
  leg->AddEntry(b4,cut3);
  leg->Draw();
  b1->GetXaxis()->SetTitle("Beta");
  /*
  int Nb=b1->GetNbinsX();
  std::cout<<"Nb = "<<Nb<<std::endl;
  TH1D *prob_beta = new TH1D("Prob","Prob",Nb,b1->GetBinLowEdge(1),b1->GetBinLowEdge(Nb)+b1->GetBinWidth(Nb));
  for(int i=1;i<=Nb;i++)
    {
      double nom =  b2->GetBinContent(i);
      double den =  b1->GetBinContent(i);
      if(den>0)prob_beta->SetBinContent(i,nom/den);
    }
  prob_beta->Draw("same");
  */
  new TCanvas("c2","c2");
  //(theta*1.74532925199432955e-02)
  t->Draw("theta>>t1",cut0,"");       
  t->Draw("theta>>t2",cut1,"same");       
  t->Draw("theta>>t3",cut2,"same");       
  t->Draw("theta>>t4",cut3,"same");       


  TH1D *t1 = (TH1D*)gROOT->FindObject("t1");
  TH1D *t2 = (TH1D*)gROOT->FindObject("t2");
  TH1D *t3 = (TH1D*)gROOT->FindObject("t3");
  TH1D *t4 = (TH1D*)gROOT->FindObject("t4");
  t1->SetLineColor(1);
  t1->SetLineStyle(2);

  t2->SetLineColor(2);
  t3->SetLineColor(3);
  t4->SetLineColor(4);

  t1->GetXaxis()->SetTitle("Theta");
  leg->Draw();

  new TCanvas("c3","c3");
  
  t->Draw("log10(peakFlux)>>pf1",cut0,"");       
  t->Draw("log10(peakFlux)>>pf2",cut1,"same");       
  t->Draw("log10(peakFlux)>>pf3",cut2,"same");       
  t->Draw("log10(peakFlux)>>pf4",cut3,"same");       


  TH1D *pf1 = (TH1D*)gROOT->FindObject("pf1");
  TH1D *pf2 = (TH1D*)gROOT->FindObject("pf2");
  TH1D *pf3 = (TH1D*)gROOT->FindObject("pf3");
  TH1D *pf4 = (TH1D*)gROOT->FindObject("pf4");
  pf1->SetLineColor(1);
  pf1->SetLineStyle(2);

  pf2->SetLineColor(2);
  pf3->SetLineColor(3);
  pf4->SetLineColor(4);

  pf1->GetXaxis()->SetTitle("log_{10}(peakFlux[50-300 keV])");
  leg->Draw();

  new TCanvas("c4","c4");
  
  t->Draw("log10(Ep)>>ep1",cut0,"");       
  t->Draw("log10(Ep)>>ep2",cut1,"same");       
  t->Draw("log10(Ep)>>ep3",cut2,"same");       
  t->Draw("log10(Ep)>>ep4","Fssc_Fsyn==0 & detected>=10","same");       


  TH1D *ep1 = (TH1D*)gROOT->FindObject("ep1");
  TH1D *ep2 = (TH1D*)gROOT->FindObject("ep2");
  TH1D *ep3 = (TH1D*)gROOT->FindObject("ep3");
  TH1D *ep4 = (TH1D*)gROOT->FindObject("ep4");
  ep1->SetLineColor(1);
  ep1->SetLineStyle(2);
  
  ep2->SetLineColor(2);
  ep3->SetLineColor(3);
  ep4->SetLineColor(4);

  ep1->GetXaxis()->SetTitle("log_{10}(Ep)");

  std::cout<<cut0.Print()<<"  "<<b1->GetEntries()<<std::endl;
  std::cout<<cut1.Print()<<" "<<b2->GetEntries()<<std::endl;
  std::cout<<cut2.Print()<<" "<<b3->GetEntries()<<std::endl;
  std::cout<<cut3.Print()<<" "<<b4->GetEntries()<<std::endl;
  int maxN=10000;
  TH1D *grbyield   = new TH1D("GRByield","GRB yield",maxN,1,maxN-1);
  TH1D *grbyield1  = new TH1D("GRByield1","GRB yield",maxN,1,maxN-1);
  double x[4]={1,10,100,1000};
  double y[4]={90,50,20,3};
  TGraph *prev  = new TGraph(4,x,y);
  grbyield->GetXaxis()->SetTitle("Number of LAT detected photons");
  grbyield->GetYaxis()->SetTitle("Number of GRB/yr");
  grbyield->SetLineColor(4);
  grbyield1->SetLineColor(2);
  grbyield->SetLineWidth(2);
  grbyield1->SetLineWidth(2);

  new TCanvas("c4","c4");
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();
  t->Draw("detected>>GRByield");
  t->Draw("detected>>GRByield1",cut0);
  std::cout<<t->GetEntries()<<" "<<GRByield->GetEntries()<<" "<<GRByield1->GetEntries()<<std::endl;

  double nrmF=t->GetEntries()/GRByield->GetEntries();
  std::cout<<nrmF<<std::endl;
  for (int i=1;i<maxN;i++)
    {
      grbyield->SetBinContent(maxN-i,nrmF*grbyield->GetBinContent(maxN-i)+grbyield->GetBinContent(maxN-i+1));
      grbyield1->SetBinContent(maxN-i,nrmF*grbyield1->GetBinContent(maxN-i)+grbyield1->GetBinContent(maxN-i+1));
    }
  grbyield->Draw("");
  grbyield1->Draw("same");
  //  grbnumber->Draw("same");
  prev->SetMarkerStyle(30);
  prev->SetMarkerColor(4);

  //  prev->Draw("*");
}
