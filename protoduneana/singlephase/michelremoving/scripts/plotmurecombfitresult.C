#include "CaloUtils.h"
#include "LanGausFit.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib/filepath_maker.h"
#include "cetlib/search_path.h"
#include "cetlib/filesystem.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TMinuit.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TLegend.h"
#include <iostream>

using namespace std;

int main(int argc, char ** argv) {

  gROOT->ProcessLine(".x /nashome/t/tjyang/.protoDUNEStyle.C");

  const int nbins = 50;
  int binsize = 10;

  bool found_fcl = false;
  string infile;
  std::string fcl_file;

  for (int iArg = 1; iArg < argc; iArg++) {
    if (!strcasecmp(argv[iArg],"-c")) {
     fcl_file = argv[++iArg];
     found_fcl = true;
    }
    if (!strcasecmp(argv[iArg],"-h")) {
      std::cout << "Usage: murecombfit " <<
                   "-c fclfile.fcl " << std::endl;
      return 1;
    }
  }

  if (!found_fcl) {
    cout << "Error: No fcl file was provided! Please provide with '-c'" << endl;
    return 0;
  }

  ////Setting up fcl parameters
  char const* fhicl_env = getenv("FHICL_FILE_PATH");
  std::string search_path;

  if (!fhicl_env) {
    std::cerr << "Expected environment variable FHICL_FILE_PATH is missing " <<
                 "or empty: using \".\"\n";
    search_path = ".";
  }
  else {
    search_path = std::string{fhicl_env};
  }

  cet::filepath_first_absolute_or_lookup_with_dot lookupPolicy{search_path};
  auto const pset = fhicl::ParameterSet::make(fcl_file, lookupPolicy);
  /////

  double alpha_orig = pset.get<double>("alpha_orig");
  double beta_orig = pset.get<double>("beta_orig");

//  double alpha_new = pset.get<double>("alpha_new");
//  double beta_new = pset.get<double>("beta_new");

  string fit_file = pset.get<string>("Fitfile");
  TFile ffit(fit_file.c_str());
  TMinuit *gMinuit = (TMinuit*)ffit.Get("gMinuit");
  double par[5];
  double epar[5];
  for (int i = 0; i<5; ++i){
    gMinuit->GetParameter(i, par[i], epar[i]);
    cout<<i<<" "<<par[i]<<" "<<epar[i]<<endl;
  }
  double alpha_new = par[3];
  double beta_new = par[4];

  std::vector<double> calconst(3);
  for (int i = 0; i<3; ++i) calconst[i] = par[i];

  string mc_file = pset.get<string>("MCfile");

  TFile fmc(mc_file.c_str());
  TTree *mctree = (TTree*)fmc.Get("calotree");
  short plane;
  float dQdx;
  float dEdx;
  float KE;
  float E;
  mctree->SetBranchAddress("plane",&plane);
  mctree->SetBranchAddress("dQdx",&dQdx);
  mctree->SetBranchAddress("dEdx",&dEdx);
  mctree->SetBranchAddress("KE",&KE);
  mctree->SetBranchAddress("E",&E);

  TH1F *dedx_mc[3][nbins];

  for (int i = 0; i<3; ++i){
    for (int j = 0; j<nbins; ++j){
      dedx_mc[i][j] = new TH1F(Form("dedx_mc_%d_%d",i,j),
                               Form("Plane %d, j %d",i,j),
                               300, 0, 15);
      dedx_mc[i][j]->Sumw2();
    }
  }

  int nentries = mctree->GetEntries();
  for (int i = 0; i<nentries; ++i){
    mctree->GetEntry(i);
    if (i%100000==0) cout<<i<<"/"<<nentries<<endl;
    int jKE = int(KE/binsize);
    if (jKE>=0 && jKE<nbins){
      dedx_mc[plane][jKE]->Fill(dEdx);
    }
  }

  double mpv[3][nbins] = {{0}};
  double errmpv[3][nbins] = {{0}};

  for (int i = 0; i<3; ++i){
    for (int j = 1; j<nbins; ++j){
      TF1 *fitsnr = runlangaufit(dedx_mc[i][j], i);
      mpv[i][j] = fitsnr->GetParameter(1);
      errmpv[i][j] = fitsnr->GetParError(1);
    }
  }       

  string data_file = pset.get<string>("Datafile");

  TFile fdata(data_file.c_str());
  TTree *datatree = (TTree*)fdata.Get("calotree");
  datatree->SetBranchAddress("plane",&plane);
  datatree->SetBranchAddress("dQdx",&dQdx);
  datatree->SetBranchAddress("dEdx",&dEdx);
  datatree->SetBranchAddress("KE",&KE);
  datatree->SetBranchAddress("E",&E);

  TH1F *dedx_data_orig[3][nbins];
  TH1F *dedx_data_new[3][nbins];

  for (int i = 0; i<3; ++i){
    for (int j = 0; j<nbins; ++j){
      dedx_data_orig[i][j] = new TH1F(Form("dedx_data_orig_%d_%d",i,j),
                                      Form("dedx_data_orig_%d_%d",i,j),
                                      300, 0, 15);
      dedx_data_new[i][j] = new TH1F(Form("dedx_data_new_%d_%d",i,j),
                                     Form("dedx_data_new_%d_%d",i,j),
                                     300, 0, 15);
    }
  }

  nentries = datatree->GetEntries();
  for (int i = 0; i<nentries; ++i){
    datatree->GetEntry(i);
    //cout<<plane<<" "<<dQdx<<" "<<dEdx<<" "<<KE<<" "<<E<<endl;
    if (i%100000==0) cout<<i<<"/"<<nentries<<endl;
    int k = int(KE/binsize);
    if (k>=0 && k<nbins){
      double dedx_orig = GetdEdx(dQdx, E, calconst[plane], alpha_orig, beta_orig);
      double dedx_new = GetdEdx(dQdx, E, calconst[plane], alpha_new, beta_new);
      dedx_data_orig[plane][k]->Fill(dedx_orig);
      dedx_data_new[plane][k]->Fill(dedx_new);
    }
  }

  double mpv_orig[3][nbins] = {};
  double errmpv_orig[3][nbins] = {};
  double mpv_new[3][nbins] = {};
  double errmpv_new[3][nbins] = {};
  double ratio_orig[3][nbins] = {};
  double ratio_new[3][nbins] = {};
  double errratio_orig[3][nbins] = {};
  double errratio_new[3][nbins] = {};
  double chi2_orig[3] = {};
  double chi2_new[3] = {};
  for (int i = 0; i<3; ++i){
    for (int j = 1; j<nbins; ++j){
      TF1 *fitsnr = runlangaufit(dedx_data_orig[i][j], i);
      mpv_orig[i][j] = fitsnr->GetParameter(1);
      errmpv_orig[i][j] = fitsnr->GetParError(1);
      fitsnr = runlangaufit(dedx_data_new[i][j], i);
      mpv_new[i][j] = fitsnr->GetParameter(1);
      errmpv_new[i][j] = fitsnr->GetParError(1);
      ratio_orig[i][j] = mpv_orig[i][j]/mpv[i][j];
      errratio_orig[i][j] = sqrt(pow(errmpv[i][j]/mpv[i][j],2)+
                                 pow(errmpv_orig[i][j]/mpv_orig[i][j],2))*ratio_orig[i][j];
      ratio_new[i][j] = mpv_new[i][j]/mpv[i][j];
      errratio_new[i][j] = sqrt(pow(errmpv[i][j]/mpv[i][j],2)+
                                pow(errmpv_new[i][j]/mpv_new[i][j],2))*ratio_new[i][j];
      chi2_orig[i] += pow(mpv_orig[i][j] - mpv[i][j],2)/(pow(errmpv[i][j],2)+pow(errmpv_orig[i][j],2));
      chi2_new[i] += pow(mpv_new[i][j] - mpv[i][j],2)/(pow(errmpv[i][j],2)+pow(errmpv_new[i][j],2));
    }
  }       

  double avgke[nbins] = {};
  for (int i = 0; i<nbins; ++i) avgke[i] = (i+0.5)*binsize;

  TGraphErrors *grmc[3], *grdata_orig[3], *grdata_new[3], *grratio_orig[3], *grratio_new[3];
  for (int i = 0; i<3; ++i){
    grmc[i] = new TGraphErrors(nbins-1, &avgke[1], &mpv[i][1], 0, &errmpv[i][1]);
    grdata_orig[i] = new TGraphErrors(nbins-1, &avgke[1], &mpv_orig[i][1], 0, &errmpv_orig[i][1]);
    grdata_new[i] = new TGraphErrors(nbins-1, &avgke[1], &mpv_new[i][1], 0, &errmpv_new[i][1]);
    grratio_orig[i] = new TGraphErrors(nbins-1, &avgke[1], &ratio_orig[i][1], 0, &errratio_orig[i][1]);
    grratio_new[i] = new TGraphErrors(nbins-1, &avgke[1], &ratio_new[i][1], 0, &errratio_new[i][1]);
  }

  TH2D *frame = new TH2D("frame",";KE (MeV);Data/MC",1000,0,500,1000,0.92,1.12);
  frame->SetStats(0);

  TCanvas *c1 = new TCanvas("c1","c1",800,800);
  TPad *pad1 = new TPad("pad1","pad1",0,0.2,1,1);
  pad1->SetBottomMargin(0.02);
  pad1->Draw();
  c1->cd();
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.2);
  pad2->SetTopMargin(0.05);
  pad2->SetBottomMargin(0.35);
  pad2->SetGridx();
  pad2->SetGridy();
  pad2->Draw();

  for (int i = 0; i<3; ++i){
    pad1->cd();
    grmc[i]->SetMarkerStyle(24);
    grmc[i]->Draw("ap");
    grmc[i]->GetXaxis()->SetRangeUser(0,500);
    grmc[i]->GetYaxis()->SetRangeUser(1,6);
    grmc[i]->GetYaxis()->SetTitle("MPV dE/dx (MeV/cm)");
    grmc[i]->SetTitle(Form("Plane %d", i));
    grdata_orig[i]->SetMarkerStyle(24);
    grdata_orig[i]->SetMarkerColor(2);
    grdata_orig[i]->SetLineColor(2);
    grdata_orig[i]->Draw("p");
    grdata_new[i]->SetMarkerStyle(24);
    grdata_new[i]->SetMarkerColor(3);
    grdata_new[i]->SetLineColor(3);
    grdata_new[i]->Draw("p");
    TPaveText *pt = new TPaveText(0.25,0.7,0.85,0.9,"NB NDC");
    pt->AddText(Form("Calibration constant = %.3f#times10^{-3}", calconst[i]*1000));
    pt->AddText(Form("Default #alpha = %.3f, #beta' = %.3f, #chi^{2} = %.0f", alpha_orig, beta_orig, chi2_orig[i]));
    pt->AddText(Form("Best fit #alpha = %.3f, #beta' = %.3f, #chi^{2} = %.0f", alpha_new, beta_new, chi2_new[i]));
    pt->Draw();
    TLegend *leg = new TLegend(0.25, 0.5, 0.85, 0.7);
    leg->SetFillStyle(0);
    leg->AddEntry(grmc[i], "MC", "p");
    leg->AddEntry(grdata_orig[i], "Data default #alpha #beta'", "p");
    leg->AddEntry(grdata_new[i], "Data best fit #alpha #beta'", "p");
    leg->Draw();

    pad2->cd();
    frame->Draw();
    frame->GetXaxis()->SetLabelSize(0.15);
    frame->GetXaxis()->SetTitleSize(0.15);
    frame->GetXaxis()->SetTitleOffset(1);
    frame->GetYaxis()->SetLabelSize(0.1);
    frame->GetYaxis()->SetTitleSize(0.15);
    frame->GetYaxis()->SetTitleOffset(.3);
    grratio_orig[i]->SetMarkerStyle(24);
    grratio_orig[i]->SetMarkerColor(2);
    grratio_orig[i]->SetLineColor(2);
    grratio_orig[i]->Draw("p");
    grratio_new[i]->SetMarkerStyle(24);
    grratio_new[i]->SetMarkerColor(3);
    grratio_new[i]->SetLineColor(3);
    grratio_new[i]->Draw("p");
    c1->Print(Form("murecombfitresult%d.png",i));
    c1->Print(Form("murecombfitresult%d.pdf",i));
  }
}
