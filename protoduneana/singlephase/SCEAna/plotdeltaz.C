#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"
#include "fhiclcpp/ParameterSet.h"
#include <string>
#include <fstream>
#include <iostream>

using namespace std;

int main(int argc, char *argv[]){

  bool found_fcl = false;
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

  string infile = pset.get<string>("infile");
  string outfile = pset.get<string>("outfile");

  float x0 = pset.get<float>("x0");
  float x1 = pset.get<float>("x1"); 

  float y0 = pset.get<float>("y0");
  float y1 = pset.get<float>("y1");
  
  TChain* shtree = new TChain("Event");
  
  std::ifstream in;
  in.open(infile.c_str());
  char line[1024];
  
  while(1){
    in.getline(line,1024);
    if (!in.good()) break;
    shtree->Add(Form("%s/michelremoving/Event", line));
  }
  in.close();
  in.clear();

  int cross_trks;
  float trkstartx[30];   //[cross_trks]
  float trkstarty[30];   //[cross_trks]
  float trkstartz[30];   //[cross_trks]
  float trkendx[30];   //[cross_trks]
  float trkendy[30];   //[cross_trks]
  float trkendz[30];   //[cross_trks]

  shtree->SetBranchAddress("cross_trks", &cross_trks);
  shtree->SetBranchAddress("trkstartx", &trkstartx);
  shtree->SetBranchAddress("trkstarty", &trkstarty);
  shtree->SetBranchAddress("trkstartz", &trkstartz);
  shtree->SetBranchAddress("trkendx", &trkendx);
  shtree->SetBranchAddress("trkendy", &trkendy);
  shtree->SetBranchAddress("trkendz", &trkendz);

  TFile *file = new TFile(outfile.c_str(), "recreate");

  TH1D *hdeltaz = new TH1D("hdeltaz", "hdeltaz", 100,0,100);
  hdeltaz->Sumw2();

  int nentries = shtree->GetEntries();
  
  for (int i = 0; i<nentries; ++i){
    if (i%10000==0) cout<<i<<"/"<<nentries<<endl;
    shtree->GetEntry(i);
    for (int j = 0; j<cross_trks; ++j){
      if (trkstartx[j]>x0 && trkstartx[j]<x1 &&
          trkstarty[j]>y0 && trkstarty[j]<y1){
        hdeltaz->Fill(trkstartz[j]);
      }
      if (trkendx[j]>x0 && trkendx[j]<x1 &&
          trkendy[j]>y0 && trkendy[j]<y1){
        hdeltaz->Fill(trkendz[j]);
      }
    }
  }
  file->Write();
  file->Close();
}
  
