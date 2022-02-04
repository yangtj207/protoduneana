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
#include "TVectorD.h"
#include "TMinuit.h"
#include "TCanvas.h"
#include "TGraph.h"
#include <iostream>

using namespace std;

const int nbins = 50;
int binsize = 10;

double pitch[nbins];

vector<float> vdedx;
vector<float> vrr;

void fcn(int& /*npar*/, double* /*gin*/, double &f, double *par, int /*iflag*/){

  double chisq = 0;
  int usedbins = 0;

  TH1F *dedx[nbins];
  for (int i = 0; i<nbins; ++i){
    dedx[i] = new TH1F(Form("dedx_%d",i),
                       Form("dedx_%d",i),
                       300, 0, 15);
  }
  double avgKE[nbins] = {};
  int nhits[nbins] = {};
  for (size_t i = 0; i<vdedx.size(); ++i){
    float ke = GetMuKEfromRange(vrr[i]+par[0]);
    int j = int(ke/binsize);
    if (j>=0 && j<nbins){
      dedx[j]->Fill(vdedx[i]);
      avgKE[j]+=ke;
      ++nhits[j];
    }
  }
  for (int i = 1; i<nbins; ++i){
    Double_t fr[2];
    Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];
    fr[0]=1.1;//this was originally 0.
    fr[1]=10.;
    if(i==1){
      fr[0]=2.2;
      fr[1]=15;
    }
    sv[0]=0.1; sv[1]=0.8*dedx[i]->GetMean(); sv[2]=dedx[i]->GetEntries()*0.05; sv[3]=0.05;
    for(int j=0; j<4; ++j){
      pllo[j] = 0.01*sv[j];
      plhi[j] = 100*sv[j];
    }
    Double_t chisqr;
    Int_t    ndf;
    Int_t    status;
    //cout<<i<<" "<<dedx[i]->GetEntries()<<endl;
    if (dedx[i]->GetEntries()<10) continue;
    TF1 *fitsnr = langaufit(dedx[i],fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf,&status);
    //cout <<"************ Fit status (FitPtr): " << status << " *********"<<endl;
    //fitsnr->SetLineColor(kRed);
    //std::cout << "************** MPV : " << fitsnr->GetParameter(1) << " +/- " << fitsnr->GetParError(1) << std::endl;
    //std::cout << "************** Chi^2/NDF : " << fitsnr->GetChisquare()/fitsnr->GetNDF() << std::endl;
    //std::cout << "   MPV : " << fitsnr->GetParameter(1) << " $$$$$$$$$$$$" << std::endl;      
    double mpv = dpdx(avgKE[i]/nhits[i], pitch[i], Mmu);
    chisq += pow(fitsnr->GetParameter(1)-mpv,2)/pow(fitsnr->GetParError(1),2);
    ++usedbins;
  }
  for (int i = 0; i<nbins; ++i){
    delete dedx[i];
  }

  f = chisq;

  cout<<"Chi2 = "<<chisq<<" "<<par[0]<<" nbins = "<<usedbins<<endl;
}
  

int main(int argc, char ** argv) {


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

  string mc_file = pset.get<string>("MCfile");

  TFile fmc(mc_file.c_str());
  TTree *mctree = (TTree*)fmc.Get("calotree");
  short plane;
  float dQdx;
  float dEdx;
  float KE;
  float RR;
  float E;
  mctree->SetBranchAddress("plane",&plane);
  mctree->SetBranchAddress("dQdx",&dQdx);
  mctree->SetBranchAddress("dEdx",&dEdx);
  mctree->SetBranchAddress("KE",&KE);
  mctree->SetBranchAddress("RR",&RR);
  mctree->SetBranchAddress("E",&E);

//  string outputfile = pset.get<string>("Outputfile");
//  TFile *foutput = new TFile(outputfile.c_str(), "recreate");

  short thisplane = pset.get<short>("thisplane");

  TVectorD *meanPitch = (TVectorD*)fmc.Get(Form("meanPitch%d",thisplane));

  for (int i = 0; i<nbins; ++i){
    pitch[i] = (*meanPitch)[i];
  }

  int nentries = mctree->GetEntries();
  for (int i = 0; i<nentries; ++i){
    mctree->GetEntry(i);
    //cout<<plane<<" "<<dQdx<<" "<<dEdx<<" "<<KE<<" "<<E<<endl;
    if (i%100000==0) cout<<i<<"/"<<nentries<<endl;
    if (plane != thisplane) continue;
    vdedx.push_back(dEdx);
    vrr.push_back(RR);
  }

  //gDirectory->cd("Rint:/");
  TMinuit *gMinuit = new TMinuit(1);
  gMinuit->SetFCN(fcn);
  double arglist[10];
  int ierflg = 0;

  // set PRINT / NO PRINT (arglist[0] = -1/0/1/2/3
  arglist[0] = 3;
  gMinuit->mnexcm("SET PRI", arglist, 1, ierflg);

  arglist[0] = 1;
  gMinuit->mnexcm("SET ERR",arglist,1,ierflg);
  arglist[0] = -1;

  /*
  double vstart = 0.4;
  double step = 0.1;
  double min = 0;
  double max = 0;
  char name[100] = "offset";
  */

  const int npt = 100;
  double vchi2[npt];
  double vshift[npt];
  for (int i = 0; i<npt; ++i){
//    vstart = i*0.1;
//    cout<<"vstart = "<<vstart<<endl;
    //gMinuit->mnparm(0,name,vstart,step,min,max,ierflg);
    double par[1] = {i*1.0/npt};
    double grad[1] = {};
    double fval;
    gMinuit->Eval(0, grad, fval, par, 1);
    cout<<par[0]<<" "<<fval<<endl;
    vchi2[i] = fval;
    vshift[i] = par[0];
    //gMinuit->mnexcm("SHOw FCNvalue", arglist, 0, ierflg);
  }

  TCanvas *c1 = new TCanvas("c1","c1");
  TGraph *gr = new TGraph(npt,vshift, vchi2);
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(0.5);
  gr->Draw("ap");
  c1->Print(Form("rrshift%d.png",thisplane));
//  arglist[0] = 500;
//  arglist[1] = 1;
//  gMinuit->mnexcm("MIGRAD",arglist,2,ierflg);

//  double par;
//  double epar;
//  gMinuit->GetParameter(0, par, epar);
//  cout<<name<<" "<<par<<" "<<epar<<endl;

  return 0;
} // main
