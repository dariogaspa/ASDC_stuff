#include "TCut.h"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFrame.h"
#include <stdio.h>
#include <string>
#include <iostream.h>
#include <fstream.h>
#include "TMath.h"
#include "TFriendElement.h"
#include "TROOT.h"


Double_t gaus1(Double_t *x, Double_t *par){
  return 1-par[0]*(1/TMath::Exp(par[1]*(x[0])));
    }





void fit_rel_nvss_30(){
  gStyle->SetOptFit( 11111);
  gStyle->SetOptStat( 00111 );
  gStyle->SetStatW( 0.2 );
  gStyle->SetStatH( 0.2 );
  gStyle->SetTitleH( 0.08 ); // Fraction of title box hight
  gStyle->SetTitleW( 0.70 ); // Fraction of title box width
  gStyle->SetTitleX( 0.1 ); // X position of the title box from left
  gStyle->SetTitleY( 1.0 ); // Y position of the title box from bottom
  gStyle->SetTitleXOffset ( 0.6 );
  gStyle->SetTitleYOffset ( 0.6 );
//TH1F*hist3 = new TH1F("hist3", "PSF R 1" ,20,-1.5, 3.05);
  Double_t dato;
  Double_t  par0gaus1_f[3];
  Int_t n =41;
  Double_t x[200];
  Double_t y[200];
  Double_t xerr[200];
  Double_t yerr[200];
  FILE  *fp;
  Double_t a[200][2];
  fp = fopen("reliability.qdp", "r");
  if (fp == 0) {
     printf("Ooops... il file non c'Ã¨\n");
     exit(1);
   }
   for(Int_t i =0;i<n; i++){
     for(Int_t k =0;k<2; k++){
       fscanf(fp, "%lf",&dato);
       a[i][k]=dato;
     }
     x[i]=a[i][0];
     //xerr[i]=a[i][1];
     y[i]=a[i][1];
     //yerr[i]=a[i][3];     
// hist3->Fill(x[i]);
     
   }


   
    
   TF1 *gaus1_f = new TF1("gaus1_f",gaus1,2,7,2);
   gaus1_f->SetParameters(par0gaus1_f[0],par0gaus1_f[1]);

   gaus1_f->SetParLimits(1,-10,10); 
   //gaus1_f->SetParLimits(4,0.5,2); 
   //   gaus1_f->FixParameter(3,1); 
   //  gaus1_f->SetParLimits(1,-1,10); 
   //   gaus1_f->FixParameter(2,2.6);
   // gaus1_f->FixParameter(3,3.2);
   // gaus1_f->FixParameter(4,0.68);   
    // gaus1_f->FixParameter(4,2.); 
   gaus1_f->SetParLimits(0,-1,1000000); 
   
   TGraph *gr1 = new TGraph(n,x,y);
 
   gr1->Fit("gaus1_f","R+","",-1.5,8);

   TCanvas *c1 = new TCanvas("c1","lightcurve",200,10,600,400);
   c1->Range(0,0,25,28);
   c1->SetFillColor( 10 );
   c1->Divide( 1, 1);
   //TGraph *gr2 = new TGraph(n1,x1,y1);
   // TCanvas *c2 = new TCanvas("c2","Two Graphs",200,10,600,400);
   // draw the graph with axis, contineous line, and put a * at each point
   c1->cd(1);
   gr1->SetMarkerColor(4);
   gr1->SetMarkerSize(1);
   gr1->SetMarkerStyle(20);

   gr1->SetLineColor(4);
   gr1->Draw("AP");
   gr1->GetXaxis()->SetTitle("LOG(LIKELIHOOD RATIO)");
   gr1->GetYaxis()->SetTitle("RELIABILITY");
   gr1->SetTitle("30GHz NVSS SAMPLE");
   double chi = gaus1_f->GetChisquare();
   cout<< chi<<endl;
   double g;
   g=0;
   for(int yf = 0; yf < 2000; yf++) {
     g=g+0.01;
     funnew=gaus1_f->Eval(g);
     //   cout<<g<<" 0 "<<funnew[g]<<" 0  "<<endl;
   }
   //gr1->

   // Hist3->fill(x);
   //  hist3->Draw();
   // gaus1_f->Draw("same");
    // superimpose the second graph by leaving out the axis option "A"
   

}
    
