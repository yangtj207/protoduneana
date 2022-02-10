#include "CaloUtils.h"
#include "LanGausFit.h"
#include "TF1.h"
#include "TSystem.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TROOT.h"
#include "TFitResult.h"
#include "TH2D.h"
#include "TVectorD.h"
#include "TGraphErrors.h"
#include "TBox.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TStyle.h"
#include <string> 
#include <fstream>

using namespace std;

/*
const int Z=18; //Atomic number of Argon
const double A=39.948; // g/mol Atomic mass of Argon
const double I=188.0e-6; // ev
const double K=0.307; // Mev.cm^2 / mol
const double Mmu=105.658; // Mev for Mu
const double Me=0.51; // Mev for electron
const double rho=1.396;//g/cm^3

double beta(double gamma){
  double value=TMath::Sqrt(1-(1.0/(gamma*gamma)));
  return value;
}

double gamma(double KE,double mass){
  double value=(double(KE)/mass)+1;
  return value;
}

double KE=266.;
double g=gamma(KE,Mmu);
double b=beta(g);

const double C=-5.2146;
const double X0=0.2;
const double X1=3.0;
const double a=0.19559;
const double m=3.0;
const double N=2*TMath::Log(10);

double density(double bg){//replaced x by x1
  double value;
  double x = TMath::Log10(bg);
  if(x<X0) return 0;
  if(x>X1) return N*x + C;
  value=a*(TMath::Power((X1-x),m));
  return N*x + C + value;
}

double Wmax(double KE,double mass){
  double g=gamma(KE,mass);
  double b=beta(g);
  double num=2*Me*(TMath::Power(b*g,2));
  double den=1+double(2*g*Me)/mass + TMath::Power((double(Me)/mass),2);
  double value=double(num)/den;
  return value;
}

double dpdx(double KE,double x,double mass){
  double g=gamma(KE,mass);
  double b=beta(g);
  double epsilon=(double(K)/2)*(double(Z)/A)*(double(x*rho)/(b*b));
  double A0=double(2*Me*(TMath::Power((b*g),2)))/I;
  double A1=double(epsilon)/I;
  double value=(1.0/x)*epsilon*((TMath::Log(A0)) + TMath::Log(A1) + 0.2 - b*b - density(b*g));
  return value;
}
*/

/*
Double_t langaufun(Double_t *x, Double_t *par) {
  Double_t invsq2pi = 0.398942280401;// Control constants
  //Double_t mpshift = -0.22278298;
  Double_t np = 500.0;
  Double_t sc = 5.0;// convolution extends to +-sc Gaussian sigmas 
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;

  //mpc = par[1]- mpshift * par[0];
  mpc=par[1];
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];
  step = (xupp-xlow)/np;
 
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }

  return (par[2] * step * sum * invsq2pi / par[3]);
}

 TF1 *langaufit(TH1 *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF, Int_t *Status)
// TF1 *langaufit(TH1D *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF)
{
  Int_t i;
  Char_t FunName[100];

  sprintf(FunName,"Fitfcn_%s",his->GetName());

  TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
  if (ffitold) delete ffitold;

  TF1 *ffit = new TF1(FunName,langaufun,fitrange[0],fitrange[1],4);
  ffit->SetParameters(startvalues);
  ffit->SetParNames("Width","MPV","Area","GSigma");
   
  for (i=0; i<4; i++) {
    ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
  }

//  his->Fit(FunName,"RB0");   // fit within specified range, use ParLimits, do not plot
  //his->SetStats(0);
  TFitResultPtr fitres = his->Fit(FunName,"RBOSQ"); // fit within specified range, use ParLimits, do not plot /////////////////// Initial code use the mode "RBO" (commented by VARUNA) ///////////

  ffit->GetParameters(fitparams);    // obtain fit parameters
  for (i=0; i<4; i++) {
    fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
  }
 
  ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
  NDF[0] = ffit->GetNDF();           // obtain ndf
  Status[0] = fitres->CovMatrixStatus();
 
  return (ffit);              // return fit function
}

/////////////////////////////// Function definition /////////////////////////////////////

Int_t langaupro(Double_t *params, Double_t &maxx, Double_t &FWHM) {
  Double_t p,x,fy,fxr,fxl;
  Double_t step;
  Double_t l,lold;
  Int_t i = 0;
  Int_t MAXCALLS = 10000;

  // Search for maximum

  p = params[1] - 0.1 * params[0];
  step = 0.05 * params[0];
  lold = -2.0;
  l = -1.0;

  while ( (l != lold) && (i < MAXCALLS) ) {
    i++;
    lold = l;
    x = p + step;
    l = langaufun(&x,params);
    if (l < lold)
      step = -step/10;
    p += step;
  }

  if (i == MAXCALLS)
    return (-1);

  maxx = x;
  fy = l/2;

  // Search for right x location of fy

  p = maxx + params[0];
  step = params[0];
  lold = -2.0;
  l    = -1e300;
  i    = 0;

  while ( (l != lold) && (i < MAXCALLS) ) {
    i++;
    lold = l;
    x = p + step;
    l = TMath::Abs(langaufun(&x,params) - fy);
    if (l > lold)
      step = -step/10;
    p += step;
  }

  if (i == MAXCALLS)
    return (-2);

  fxr = x;

  // Search for left x location of fy

  p = maxx - 0.5 * params[0];
  step = -params[0];
  lold = -2.0;
  l    = -1e300;
  i    = 0;

  while ( (l != lold) && (i < MAXCALLS) ) {
    i++;
    lold = l;
    x = p + step;
    l = TMath::Abs(langaufun(&x,params) - fy);
    if (l > lold)
      step = -step/10;
    p += step;
  }

  if (i == MAXCALLS)
    return (-3);

  fxl = x;
  FWHM = fxr - fxl;
  return (0);
}
*/

//void fit(){
int main(int argc, char *argv[]) {

  gROOT->ProcessLine(".x /nashome/t/tjyang/.protoDUNEStyle.C");
  if (!argv[1]){
    cout<<"Error: no input file"<<endl;
  }
  gStyle->SetOptStat(2210);
  //string filename = "Validate_mich2_r5387.root";
  string filename = argv[1];

  TFile f(filename.c_str());

  double calib_const[3];
  ifstream in;
  char line[1024];
  for (int i = 0; i<3; ++i){
    in.open(Form("calconst_%d.txt", i));
    while(1){
      in.getline(line,1024);
      if(!in.good()) break;
      double error;
      sscanf(line, "%lf %lf", &calib_const[i], &error);
    }
    in.close();
    in.clear();
  }

  const size_t nbins = 50;
  int binsize = 10;

  TH1D *dedx[3][nbins];
  TH2D *dedxke[3];
  //TH1F *hdqdx_mip[3];
  TH1F *hdedx_mip[3];

  TCanvas *can = new TCanvas("can","can");
  TCanvas *can2 = new TCanvas("can2","can2",800,800);
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.2, 1, 1.);
  pad1->SetBottomMargin(0.02);
  //  pad1->SetGridx();
  //  pad1->SetGridy();
  pad1->Draw();
  can2->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.2);
  pad2->SetTopMargin(0.05);
  pad2->SetBottomMargin(0.35);
  pad2->SetGridx();
  pad2->SetGridy();
  pad2->Draw();

  TH2D *frame = new TH2D("frame",";KE (MeV);Data/MC",1000,0,500,1000,0.85,1.05);
  frame->SetStats(0);

//  TLatex latex;
//  latex.SetTextSize(0.04);
//  latex.SetTextColor(2);

  TVectorD *params = (TVectorD*)f.Get("params");

  for (size_t i = 0; i<3; ++i){
    can->cd();
    can->Print(Form("dedx_%zu.ps[", i));
    dedxke[i] = (TH2D*)f.Get(Form("dedxke_%zu",i));
    //hdqdx_mip[i] = (TH1F*)f.Get(Form("hdqdx_mip_%zu",i));
    hdedx_mip[i] = (TH1F*)f.Get(Form("hdedx_mip_%zu",i));
    vector<double> vke;
    vector<double> vdedx;
    vector<double> ededx;
    vector<double> vdedxth;
    vector<double> ratio;
    vector<double> eratio;
    vector<double> lwidth;
    vector<double> elwidth;
    vector<double> gwidth;
    vector<double> egwidth;
    TVectorD *meanKE = (TVectorD*)f.Get(Form("meanKE%zu",i));
    TVectorD *meanPitch = (TVectorD*)f.Get(Form("meanPitch%zu",i));
    double totalchi2 = 0;
    double chi2_250_450 = 0;
    double avgpitch = 0;
    for (size_t j = 1; j<nbins; ++j){
      dedx[i][j] = (TH1D*)f.Get(Form("dedx_%zu_%zu", i, j));
      dedx[i][j]->Draw();
      /*
      Double_t fr[2];
      Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];
      fr[0]=1.1;//this was originally 0.
      fr[1]=10.;
      if(j==1){
        fr[0]=2.2;
        fr[1]=15;
      }
      if (dedx[i][j]->GetMean()<10){
        //sv[0]=0.13*dedx[i]->GetRMS(); sv[1]=0.8*dedx[i]->GetMean(); sv[2]=dedx[i]->GetEntries()*0.1; sv[3]=.3;
        if (i==2){
          sv[0]=0.08; sv[1]=0.8*dedx[i][j]->GetMean(); sv[2]=dedx[i][j]->GetEntries()*0.05; sv[3]=0.13;
        }
        else{
          sv[0]=0.08; sv[1]=0.8*dedx[i][j]->GetMean(); sv[2]=dedx[i][j]->GetEntries()*0.05; sv[3]=0.15;
        }

        //if(j==0){ sv[0]=0.2; sv[1]=4.7; sv[2]=20; sv[3]=.01;}
        //if(j==1){ sv[0]=0.2; sv[1]=3.0; sv[2]=10; sv[3]=.01;}
        //if(j==2){ sv[1]=2.5;}
        //if(j==3){ sv[1]=2.0;}
        //if(j==4){ sv[1]=2.0;}
      }
      else{
        sv[0]=0.16*dedx[i][j]->GetRMS(); sv[1]=0.9*dedx[i][j]->GetMean(); sv[2]=dedx[i][j]->GetEntries()*100; sv[3]=dedx[i][j]->GetRMS()/5.;
      }
      for(int k=0; k<4; ++k){
        pllo[k] = 0.01*sv[k];
        plhi[k] = 100*sv[k];
      }
      Double_t chisqr;
      Int_t    ndf;
      Int_t    status;
      //    TF1 *fitsnr = langaufit(dedx[i][j],fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
      */
      //TF1 *fitsnr = langaufit(dedx[i][j],fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf,&status);
      TF1 *fitsnr = runlangaufit(dedx[i][j], i);
      //    cout <<"************ Fit status (gMinuit): " << gMinuit << ", "<< gMinuit->fCstatu.Data() <<" *********"<<endl;
      //cout <<"************ Fit status (FitPtr): " << status << " *********"<<endl;
      fitsnr->SetLineColor(kRed);
      std::cout << "************** MPV : " << fitsnr->GetParameter(1) << " +/- " << fitsnr->GetParError(1) << std::endl;
      std::cout << "************** Chi^2/NDF : " << fitsnr->GetChisquare()/fitsnr->GetNDF() << std::endl;
      std::cout << "   MPV : " << fitsnr->GetParameter(1) << " $$$$$$$$$$$$" << std::endl;
      
      can->Print(Form("dedx_%zu.ps", i));

      vke.push_back((*meanKE)[j]);
      vdedx.push_back(fitsnr->GetParameter(1));
      ededx.push_back(fitsnr->GetParError(1));
      lwidth.push_back(fitsnr->GetParameter(0));
      elwidth.push_back(fitsnr->GetParError(0));
      gwidth.push_back(fitsnr->GetParameter(3));
      egwidth.push_back(fitsnr->GetParError(3));
      vdedxth.push_back(dpdx((*meanKE)[j], (*meanPitch)[j], Mmu));
      ratio.push_back(vdedx.back()/vdedxth.back());
      eratio.push_back(ededx.back()/vdedxth.back());
      totalchi2 += pow((vdedx.back()-vdedxth.back())/ededx.back(),2);
      if ((j+0.5)*binsize > 250 && (j+0.5)*binsize < 450){
        chi2_250_450 += pow((vdedx.back()-vdedxth.back())/ededx.back(),2);
      }
      avgpitch += (*meanPitch)[j];
    }
    avgpitch /= nbins - 1;
    can->Print(Form("dedx_%zu.ps]", i));

    //Fit MIP dE/dx distributions
    Double_t fr[2];
    Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];
    fr[0]=1.1;//this was originally 0.
    fr[1]=10.;
    sv[0]=0.1; sv[1]=0.8*hdedx_mip[i]->GetMean(); sv[2]=hdedx_mip[i]->GetEntries()*0.05; sv[3]=0.05;
    for(int k=0; k<4; ++k){
      pllo[k] = 0.01*sv[k];
      plhi[k] = 100*sv[k];
    }
    Double_t chisqr;
    Int_t    ndf;
    Int_t    status;
    TF1 *fitsnr = langaufit(hdedx_mip[i],fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf,&status);
    cout <<"************ Fit status (FitPtr): " << status << " *********"<<endl;
    fitsnr->SetLineColor(kRed);
    std::cout << "************** MPV : " << fitsnr->GetParameter(1) << " +/- " << fitsnr->GetParError(1) << std::endl;
    std::cout << "************** Chi^2/NDF : " << fitsnr->GetChisquare()/fitsnr->GetNDF() << std::endl;
    std::cout << "   MPV : " << fitsnr->GetParameter(1) << " $$$$$$$$$$$$" << std::endl;
    can->Print(Form("dedx_mip_%zu.pdf",i));
    can->Print(Form("dedx_mip_%zu.png",i));

    pad1->cd();
    TGraphErrors *gr = new TGraphErrors(vke.size(), &vke[0], &vdedx[0], 0, &ededx[0]);
    TGraph *grth = new TGraph(vke.size(), &vke[0], &vdedxth[0]);
    TGraphErrors *grratio = new TGraphErrors(vke.size(), &vke[0], &ratio[0], 0, &eratio[0]);
    dedxke[i]->SetStats(0);
    dedxke[i]->Draw("colz");
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(2);
    gr->SetMarkerSize(0.5);
    gr->Draw("pe");
    grth->SetLineWidth(2);
    grth->SetLineColor(2);
    grth->Draw("c");
    TPaveText *pt = new TPaveText(0.25,0.7,0.85,0.9,"NB NDC");
    pt->AddText(Form("Calibration constant = %.3f#times10^{-3}", calib_const[i]*1000));
    pt->AddText(Form("#alpha = %.3f, #beta' = %.3f", (*params)[0], (*params)[1]));
    pt->AddText(Form("Average track pitch = %.2f", avgpitch));
    pt->AddText(Form("For 50<KE<500 MeV, #chi^{2} = %.1f, nbins = %zu", totalchi2, nbins-1));
    pt->AddText(Form("For 250<KE<450 MeV, #chi^{2} = %.1f, nbins = %d", chi2_250_450, 200/binsize));
    pt->Draw();
    pad2->cd();
    frame->Draw();
    frame->GetXaxis()->SetLabelSize(0.15);
    frame->GetXaxis()->SetTitleSize(0.15);
    frame->GetXaxis()->SetTitleOffset(1);
    frame->GetYaxis()->SetLabelSize(0.1);
    frame->GetYaxis()->SetTitleSize(0.15);
    frame->GetYaxis()->SetTitleOffset(.3);
    TBox *box = new TBox(250,0.85,450,1.05);
    box->SetFillColorAlpha(kRed,0.2);
    box->SetLineColor(kRed);
    box->Draw();
    grratio->Draw("p");
    can2->Print(Form("dedxke_%zu.png",i));
    can2->Print(Form("dedxke_%zu.pdf",i));
    can->cd();
    TGraphErrors *grlwidth = new TGraphErrors(vke.size(), &vke[0], &lwidth[0], 0, &elwidth[0]);
    TGraphErrors *grgwidth = new TGraphErrors(vke.size(), &vke[0], &gwidth[0], 0, &egwidth[0]);
    grlwidth->SetMarkerStyle(20);
    grlwidth->SetTitle(Form("Plane %zu", i));
    grlwidth->GetXaxis()->SetTitle("KE (MeV)");
    grlwidth->GetYaxis()->SetTitle("Landau width");
    grlwidth->Draw("ap");
    grlwidth->Fit("pol2");
    can->Print(Form("lwidth_%zu.png",i));
    can->Print(Form("lwidth_%zu.pdf",i));
    grgwidth->SetMarkerStyle(20);
    grgwidth->SetTitle(Form("Plane %zu", i));
    grgwidth->GetXaxis()->SetTitle("KE (MeV)");
    grgwidth->GetYaxis()->SetTitle("Gaussian width");
    grgwidth->Draw("ap");
    grgwidth->Fit("pol2");
    can->Print(Form("gwidth_%zu.png",i));
    can->Print(Form("gwidth_%zu.pdf",i));
  }

}

  
