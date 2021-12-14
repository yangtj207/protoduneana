#include <string>
#include <iostream>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TLatex.h>
#include <TF1.h>
#include <TLine.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TROOT.h>
#include "TVectorD.h"
#include <fstream> 

using namespace std;

int main(int argc, char *argv[]) {

  gROOT->ProcessLine(".x /nashome/t/tjyang/.protoDUNEStyle.C");

  if (!argv[1]){
    cout<<"Error: no input file"<<endl;
  }

  string filename = argv[1];

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TCanvas *can[3];
  
  TFile f(filename.c_str());

  TVectorD *NDF = (TVectorD*)f.Get("NDF");
  TVectorD *meanPitch = (TVectorD*)f.Get("meanPitch");

  TGraph *gr[3];

  TLatex tt;
  tt.SetNDC();
  tt.SetTextSize(0.04);

  for (int i = 0; i<3; ++i){
    
    can[i] = new TCanvas(Form("can_%d",i), Form("can_%d",i));
    gr[i] = (TGraph*)f.Get(Form("chi2_plane_%d",i));
    gr[i]->SetMarkerStyle(20);
    gr[i]->SetMarkerSize(0.5);
    gr[i]->SetTitle(Form("Plane %d",i));
    gr[i]->GetXaxis()->SetTitle("Calibration constant");
    gr[i]->GetYaxis()->SetTitle("#chi^{2}");
    double bestx = 0;
    double miny = 1e10;
    for (int j = 0; j<gr[i]->GetN(); ++j){
      if (gr[i]->GetPointY(j)<miny){
        miny = gr[i]->GetPointY(j);
        bestx = gr[i]->GetPointX(j);
      }
    }
    gr[i]->Draw("ap");
    //gr[i]->GetXaxis()->SetRangeUser(bestx-0.002e-3, bestx+0.002e-3);
    //gr[i]->GetYaxis()->SetRangeUser(0, miny+10);
    gr[i]->Fit("pol2","RQ","",bestx-0.002e-3, bestx+0.002e-3);
    TF1 *fun = (TF1*)gr[i]->FindObject("pol2");
    double a = fun->GetParameter(0);
    double b = fun->GetParameter(1);
    double c = fun->GetParameter(2);
    double C0 = -b/(2*c);
    double chisq0 = fun->Eval(C0);
    double C1 = (-b-sqrt(b*b-4*c*(a-chisq0-1)))/(2*c);
    double C2 = (-b+sqrt(b*b-4*c*(a-chisq0-1)))/(2*c);
    double x1 = C0-(C0-C1)*3;
    double x2 = C0+(C0-C1)*3;
    double y1 = fun->Eval(x1);
    double y2 = fun->Eval(x2);
    gr[i]->GetXaxis()->SetRangeUser(x1,x2);
    gr[i]->GetYaxis()->SetRangeUser(miny-1,TMath::Max(y1,y2));
    TLine *l1 = new TLine(C1,miny-1,C1,fun->Eval(C1));
    l1->SetLineStyle(2);
    l1->Draw();
    TLine *l2 = new TLine(C1,fun->Eval(C1),C2,fun->Eval(C2));
    l2->SetLineStyle(2);
    l2->Draw();
    TLine *l3 = new TLine(C2,miny-1,C2,fun->Eval(C2));
    l3->SetLineStyle(2);
    l3->Draw();
    tt.DrawLatex(0.3, 0.8, Form("#chi^{2}_{min}:%.2f, nbins:%d", chisq0, int((*NDF)[i])));
    if (i!=2){
      tt.DrawLatex(0.3, 0.72, Form("Calib_const = (%.3f#pm%.3f)#times10^{-3}", C0*1e3, (C0-C1)*1e3));
    }
    else{
      tt.DrawLatex(0.3, 0.72, Form("Calib_const = (%.4f#pm%.4f)#times10^{-3}", C0*1e3, (C0-C1)*1e3));
    }
    tt.DrawLatex(0.3, 0.64, Form("Mean pitch: %.3f", (*meanPitch)[i]));
    //cout<<a<<" "<<b<<" "<<c<<" "<<chisq0<<" "<<C0<<"-"<<C0-C1<<"+"<<C2-C0<<endl;
    cout<<"Plane "<<i<<endl;
    cout<<"Minimal chi2 "<<chisq0<<endl;
    cout<<"Calibration constant "<<C0<<"+-"<<C0-C1<<endl;
    can[i]->Print(Form("calconst_%d.png",i));
    can[i]->Print(Form("calconst_%d.pdf",i));
    ofstream outfile;
    outfile.open(Form("calconst_%d.txt",i));
    outfile<<C0<<" "<<C0-C1<<endl;
    outfile.close();
  }

}
