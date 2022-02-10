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
#include <iostream>

using namespace std;

const int nbins = 50;
int binsize = 10;

vector<vector<double>> dqdx_data(3);
vector<vector<double>> E_data(3);
vector<vector<double>> ke_data(3);

double mpv[3][nbins] = {{0}};
double errmpv[3][nbins] = {{0}};

void fcn(int& /*npar*/, double* /*gin*/, double &f, double *par, int /*iflag*/){

  double chisq = 0;

  TH1D *dedx_data[3][nbins];
  for (int i = 0; i<3; ++i){
    for (int j = 0; j<nbins; ++j){
      dedx_data[i][j] = new TH1D(Form("dedx_data_%d_%d",i,j),
                                 Form("dedx_data_%d_%d",i,j),
                                 300, 0, 15);
    }
    for (size_t j = 0; j<dqdx_data[i].size(); ++j){
      int k = int(ke_data[i][j]/binsize);
      if (k>=0 && k<nbins){
        double dedx = GetdEdx(dqdx_data[i][j], E_data[i][j], par[i], par[3], par[4]);
        dedx_data[i][k]->Fill(dedx);
      }
    }
    for (int j = 1; j<nbins; ++j){
      /*
      Double_t fr[2];
      Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];
      fr[0]=1.1;//this was originally 0.
      fr[1]=10.;
      if(j==1){
        fr[0]=2.2;
        fr[1]=15;
      }
      sv[0]=0.1; sv[1]=0.8*dedx_data[i][j]->GetMean(); sv[2]=dedx_data[i][j]->GetEntries()*0.05; sv[3]=0.05;
      for(int k=0; k<4; ++k){
        pllo[k] = 0.01*sv[k];
        plhi[k] = 100*sv[k];
      }
      Double_t chisqr;
      Int_t    ndf;
      Int_t    status;
      */
      TF1 *fitsnr = runlangaufit(dedx_data[i][j], i);
      //cout <<"************ Fit status (FitPtr): " << status << " *********"<<endl;
      //fitsnr->SetLineColor(kRed);
      //std::cout << "************** MPV : " << fitsnr->GetParameter(1) << " +/- " << fitsnr->GetParError(1) << std::endl;
      //std::cout << "************** Chi^2/NDF : " << fitsnr->GetChisquare()/fitsnr->GetNDF() << std::endl;
      //std::cout << "   MPV : " << fitsnr->GetParameter(1) << " $$$$$$$$$$$$" << std::endl;      
      chisq += pow(fitsnr->GetParameter(1)-mpv[i][j],2)/(pow(errmpv[i][j],2)+pow(fitsnr->GetParError(1),2));
    }
    for (int j = 0; j<nbins; ++j){
      delete dedx_data[i][j];
    }
  }
  cout<<"Chi2 = "<<chisq<<" "<<par[0]<<" "<<par[1]<<" "<<par[2]<<" "<<par[3]<<" "<<par[4]<<endl;
  f = chisq;
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
  float E;
  mctree->SetBranchAddress("plane",&plane);
  mctree->SetBranchAddress("dQdx",&dQdx);
  mctree->SetBranchAddress("dEdx",&dEdx);
  mctree->SetBranchAddress("KE",&KE);
  mctree->SetBranchAddress("E",&E);

  string outputfile = pset.get<string>("Outputfile");
  TFile *foutput = new TFile(outputfile.c_str(), "recreate");

  //TH2D *dedxke_data[3];
  TH2D *dedxke_mc[3];

  for (int i = 0; i < 3; ++i) {
//    dedxke_data[i] = new TH2D(Form("dedxke_data_%d", i),
//                              Form("Data, plane:%d;KE (MeV);dE/dx (MeV/cm)", i),
//                              250,0,500,80,0,8);
    dedxke_mc[i] = new TH2D(Form("dedxke_mc_%d", i),
                              Form("MC, plane:%d;KE (MeV);dE/dx (MeV/cm)", i),
                              250,0,500,80,0,8);
  }

  //TH1D *dedx_data[3][nbins];
  TH1D *dedx_mc[3][nbins];

  for (int i = 0; i<3; ++i){
    for (int j = 0; j<nbins; ++j){
//      dedx_data[i][j] = new TH1D(Form("dedx_data_%d_%d",i,j),
//                                 Form("Plane %d, j %d",i,j),
//                                 300, 0, 15);
//      dedx_data[i][j]->Sumw2();
      dedx_mc[i][j] = new TH1D(Form("dedx_mc_%d_%d",i,j),
                               Form("Plane %d, j %d",i,j),
                               300, 0, 15);
      dedx_mc[i][j]->Sumw2();
    }
  }

  int nentries = mctree->GetEntries();
  for (int i = 0; i<nentries; ++i){
    mctree->GetEntry(i);
    //cout<<plane<<" "<<dQdx<<" "<<dEdx<<" "<<KE<<" "<<E<<endl;
    if (i%100000==0) cout<<i<<"/"<<nentries<<endl;
    dedxke_mc[plane]->Fill(KE, dEdx);
    int jKE = int(KE/binsize);
    if (jKE>=0 && jKE<nbins){
      dedx_mc[plane][jKE]->Fill(dEdx);
    }
  }

  for (int i = 0; i<3; ++i){
    for (int j = 1; j<nbins; ++j){
      /*
      Double_t fr[2];
      Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];
      fr[0]=1.1;//this was originally 0.
      fr[1]=10.;
      if(j==1){
        fr[0]=2.2;
        fr[1]=15;
      }
      sv[0]=0.1; sv[1]=0.8*dedx_mc[i][j]->GetMean(); sv[2]=dedx_mc[i][j]->GetEntries()*0.05; sv[3]=0.05;
      for(int k=0; k<4; ++k){
        pllo[k] = 0.01*sv[k];
        plhi[k] = 100*sv[k];
      }
      Double_t chisqr;
      Int_t    ndf;
      Int_t    status;
      */
      TF1 *fitsnr = runlangaufit(dedx_mc[i][j], i);
      //cout <<"************ Fit status (FitPtr): " << status << " *********"<<endl;
      //fitsnr->SetLineColor(kRed);
//      std::cout << "************** MPV : " << fitsnr->GetParameter(1) << " +/- " << fitsnr->GetParError(1) << std::endl;
//      std::cout << "************** Chi^2/NDF : " << fitsnr->GetChisquare()/fitsnr->GetNDF() << std::endl;
//      std::cout << "   MPV : " << fitsnr->GetParameter(1) << " $$$$$$$$$$$$" << std::endl;
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

  nentries = datatree->GetEntries();
  for (int i = 0; i<nentries; ++i){
    datatree->GetEntry(i);
    //cout<<plane<<" "<<dQdx<<" "<<dEdx<<" "<<KE<<" "<<E<<endl;
    if (i%100000==0) cout<<i<<"/"<<nentries<<endl;
    dqdx_data[plane].push_back(dQdx);
    E_data[plane].push_back(E);
    ke_data[plane].push_back(KE);
  }

  //gDirectory->cd("Rint:/");
  TMinuit *gMinuit = new TMinuit(5);
  gMinuit->SetFCN(fcn);
  double arglist[10];
  int ierflg = 0;

  arglist[0] = 1;
  gMinuit->mnexcm("SET ERR",arglist,1,ierflg);

  double vstart[5] = {1.166e-3,
                      1.122e-3,
                      1.038e-3,
                      0.93,
                      0.212};
  double step[5] = {0.1e-3,
                    0.1e-3,
                    0.1e-3,
                    0.01,
                    0.01};
  double min[5] = {0.9e-3,
                   0.9e-3,
                   0.9e-3,
                   0.,
                   0.};
  double max[5] = {1.4e-3,
                   1.4e-3,
                   1.4e-3,
                   0,
                   0.};
  char name[5][100] = {"calconst0",
                       "calconst1",
                       "calconst2",
                       "alpha",
                       "beta"};
  
  for (int i = 0; i<5; ++i){
    gMinuit->mnparm(i,name[i],vstart[i],step[i],min[i],max[i],ierflg);
  }
  gMinuit->FixParameter(0);
  gMinuit->FixParameter(1);
  gMinuit->FixParameter(2);

  gMinuit->mnexcm("show eps", arglist, 0, ierflg);

  //set EPS precision
  arglist[0] = 1.e-7;
  gMinuit->mnexcm("SET EPS",arglist,1,ierflg);

  //set strategy
//  arglist[0] = 0;
//  gMinuit->mnexcm("SET STR",arglist,1,ierflg);

  arglist[0] = 5;
  arglist[1] = 10;
  arglist[2] = 0.2;
  arglist[3] = 0.25;
  gMinuit->mnexcm("scan", arglist, 4, ierflg);

  arglist[0] = 4;
  arglist[1] = 10;
  arglist[2] = 0.925;
  arglist[3] = 0.935;
  gMinuit->mnexcm("scan", arglist, 4, ierflg);

  arglist[0] = 500;
  arglist[1] = 1;
  //gMinuit->mnexcm("MIGRAD",arglist,2,ierflg);
  gMinuit->mnexcm("MINIMIZE",arglist,2,ierflg);

  double par[5];
  double epar[5];
  for (int i = 0; i<5; ++i){
    gMinuit->GetParameter(i, par[i], epar[i]);
    cout<<name[i]<<" "<<par[i]<<" "<<epar[i]<<endl;
  }

  foutput->cd();
  foutput->Write();
  gMinuit->Write("gMinuit");
  foutput->Close();

  return 0;
} // main
