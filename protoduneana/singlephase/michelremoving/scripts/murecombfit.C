#include "LanGausFit.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib/filepath_maker.h"
#include "cetlib/search_path.h"
#include "cetlib/filesystem.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include <iostream>

using namespace std;

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

  const int nbins = 50;
  int binsize = 10;

  //TH1F *dedx_data[3][nbins];
  TH1F *dedx_mc[3][nbins];

  for (int i = 0; i<3; ++i){
    for (int j = 0; j<nbins; ++j){
//      dedx_data[i][j] = new TH1F(Form("dedx_data_%d_%d",i,j),
//                                 Form("Plane %d, j %d",i,j),
//                                 300, 0, 15);
//      dedx_data[i][j]->Sumw2();
      dedx_mc[i][j] = new TH1F(Form("dedx_mc_%d_%d",i,j),
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
      TF1 *fitsnr = langaufit(dedx_mc[i][j],fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf,&status);
      cout <<"************ Fit status (FitPtr): " << status << " *********"<<endl;
      fitsnr->SetLineColor(kRed);
      std::cout << "************** MPV : " << fitsnr->GetParameter(1) << " +/- " << fitsnr->GetParError(1) << std::endl;
      std::cout << "************** Chi^2/NDF : " << fitsnr->GetChisquare()/fitsnr->GetNDF() << std::endl;
      std::cout << "   MPV : " << fitsnr->GetParameter(1) << " $$$$$$$$$$$$" << std::endl;
    }
  }       

  foutput->Write();
  foutput->Close();

  return 0;
} // main
