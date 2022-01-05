#include <TFile.h>
#include <TH1F.h>
#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char ** argv) {

  if (!argv[1]) {
    cout << "Error: No input file was provided!" << endl;
    cout << "Usage: " << endl;
    cout << "ConvertXCorrection root_file_or_list" << endl;

    return 0;
  }
  
  string infile = argv[1];

  int position1 = infile.find(".root");
  int position2 = infile.find("_r");
  int tv = std::stoi(infile.substr(position2+2,position1-position2-2));
  std::cout<<tv<<std::endl;
  TFile* fin = new TFile(infile.c_str());
  for (int plane = 0; plane < 3; ++plane){
    TH1F* hCorr = (TH1F*)fin->Get(Form("dqdx_X_correction_hist_%i",plane));
    
    //  hCorr->Draw();
    
    int nbins = hCorr->GetXaxis()->GetNbins();
    float x; 
    float dx;
    float corr;
    std::ofstream fout;
    fout.open(Form("xcorr_%d_%d.csv",tv,plane));
    //int tv = 0;
    
    //fout << "tv,channel,x,dx,shape,shape_err" << std::endl;
    fout << "channel,tv,x,dx,shape,shape_err" << std::endl;
    
    for (int i=1; i<=nbins; ++i) {
      x = hCorr->GetBinCenter(i);
      dx = hCorr->GetBinWidth(i);
      corr = hCorr->GetBinContent(i);
      if (corr == 0.) corr = 1.0;
      //fout << tv << "," << plane*10000+i << "," << x << "," << dx << "," << corr << "," << corr*0.1 << std::endl;
      fout << plane*10000+i << "," << tv << "," << x << "," << dx << "," << corr << "," << corr*0.1 << std::endl;
    }
    
    fout.close();
  }
}
