#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char ** argv) {

  if (!argv[1]) {
    cout << "Error: No input file was provided!" << endl;
    cout << "Usage: " << endl;
    cout << "ConvertYZCorrection root_file_or_list" << endl;

    return 0;
  }

  string infile = argv[1];

  int position1 = infile.find(".root");
  int position2 = infile.find("_r");
  int tv = std::stoi(infile.substr(position2+2,position1-position2-2));
  std::cout<<tv<<std::endl;
  TFile* fin = new TFile(infile.c_str());
  for (int plane = 0; plane < 3; ++plane){
    TH2F* hCorr[2];
    hCorr[0] = (TH2F*)fin->Get(Form("correction_dqdx_ZvsY_negativeX_hist_%i",plane));
    hCorr[1] = (TH2F*)fin->Get(Form("correction_dqdx_ZvsY_positiveX_hist_%i",plane));
    
    std::ofstream fout;
    fout.open(Form("yzcorr_%d_%d.csv",tv,plane));
    //int tv = 0;
    float y, z; 
    float dy, dz;
    float corr;
    
    //fout << "tv,channel,y,dy,z,dz,corr,corr_err" << std::endl;
    fout << "channel,tv,y,dy,z,dz,corr,corr_err" << std::endl;
    for (int side=0; side<2; ++side) {
      int nbinsX = hCorr[side]->GetXaxis()->GetNbins();
      int nbinsY = hCorr[side]->GetYaxis()->GetNbins();
      
      for (int i=1; i<=nbinsX; ++i) {
        z = hCorr[side]->GetXaxis()->GetBinCenter(i);
        dz = hCorr[side]->GetXaxis()->GetBinWidth(i);
        for (int j=1; j<=nbinsY; ++j) {
          y = hCorr[side]->GetYaxis()->GetBinCenter(j);
          dy = hCorr[side]->GetYaxis()->GetBinWidth(j);
          corr = hCorr[side]->GetBinContent(i,j);
          if (corr == 0.) corr = 1.0;
          //fout << tv << "," << plane*10000000+side*1000000+i*1000+j << "," << y << "," << dy << "," << z << "," << dz << "," << corr << "," << corr*0.1 << std::endl;
          fout << plane*10000000+side*1000000+i*1000+j << "," << tv << "," << y << "," << dy << "," << z << "," << dz << "," << corr << "," << corr*0.1 << std::endl;
        }
      }
    }

    fout.close();
  }
}
