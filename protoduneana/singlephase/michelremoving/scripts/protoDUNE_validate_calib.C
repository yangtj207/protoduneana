#define protoDUNE_validate_calib_cxx
#include "protoDUNE_validate_calib.h"
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TStyle.h>
#include <iostream>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TMath.h>
#include <TSpline.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TVector3.h>
#include <TVectorD.h>
#include <fstream>
#include <vector>
#include <cmath>
#include <math.h> 
#include <string>
#include <TImage.h>
#include <iomanip>
#include <algorithm>

#include "fhiclcpp/ParameterSet.h"

using namespace std;

////defining recombination function
float LAr_density=1.39;
//ArgoNeuT parameters
float alp=0.93;
float bet=0.212;
//Abbey's parameters https://indico.fnal.gov/event/46502/contributions/206721/attachments/139268/174710/recombination_20210126.pdf
//float alp=0.854;
//float bet=0.208;
//Abbey's parameters https://indico.fnal.gov/event/46503/contributions/215375/attachments/143230/181102/recombination_20210519.pdf
//float alp=0.912;
//float bet=0.195;
//float dedx=2.08;
//bool userecom=true;
bool recalib=true;
bool sceon=true;
bool corr_end;
int xbins;
double xmin;
double xmax;
vector<double> rroffset;
string outfile;
//float recom_factor(float totEf){
//  if (!userecom) return 1;
//  float xsi=bet*dedx/(LAr_density*totEf);
//  float xsi0=bet*dedx/(LAr_density*0.4867);
//  float rec0=log(alp+xsi0)/xsi0;
//  return (rec0*xsi)/log(alp+xsi);
//}

float GetdEdx(float dQdx, float E_field){
  float Beta = bet/(LAr_density*E_field);
  double Wion = 23.6e-6;
  return (exp(Beta * Wion *dQdx) - alp) / Beta;
}

TFile *ef = TFile::Open("$DUNE_PARDATA_DIR/SpaceChargeProtoDUNE/SCE_DataDriven_180kV_v4.root");
TH3F *xneg=(TH3F*)ef->Get("Reco_ElecField_X_Neg");
TH3F *yneg=(TH3F*)ef->Get("Reco_ElecField_Y_Neg");
TH3F *zneg=(TH3F*)ef->Get("Reco_ElecField_Z_Neg");
TH3F *xpos=(TH3F*)ef->Get("Reco_ElecField_X_Pos");
TH3F *ypos=(TH3F*)ef->Get("Reco_ElecField_Y_Pos");
TH3F *zpos=(TH3F*)ef->Get("Reco_ElecField_Z_Pos");
float tot_Ef(float xval,float yval,float zval){
  float E0value=0.4867;
  if (!sceon){
    return E0value;
  }
  else{
    if(xval>=0){
      float ex=E0value+E0value*xpos->GetBinContent(xpos->FindBin(xval,yval,zval));
      float ey=0.0+E0value*ypos->GetBinContent(ypos->FindBin(xval,yval,zval));
      float ez=0.0+E0value*zpos->GetBinContent(zpos->FindBin(xval,yval,zval));
      return sqrt(ex*ex+ey*ey+ez*ez);
      // return ex;
    }
    else{
      float ex=E0value+E0value*xneg->GetBinContent(xneg->FindBin(xval,yval,zval));
      float ey=0.0+E0value*yneg->GetBinContent(yneg->FindBin(xval,yval,zval));
      float ez=0.0+E0value*zneg->GetBinContent(zneg->FindBin(xval,yval,zval));
      return sqrt(ex*ex+ey*ey+ez*ez);
      // return ex;
    }
  }
}

TH3F *zpos_fd=(TH3F*)ef->Get("RecoFwd_Displacement_Z_Pos");
TH3F *zneg_fd=(TH3F*)ef->Get("RecoFwd_Displacement_Z_Neg");
TH3F *zpos_bd=(TH3F*)ef->Get("RecoBkwd_Displacement_Z_Pos");
TH3F *zneg_bd=(TH3F*)ef->Get("RecoBkwd_Displacement_Z_Neg");

float zoffset(float xval,float yval,float zval){
  if(xval>=0){
    return zpos_fd->GetBinContent(zpos_fd->FindBin(xval,yval,zval));
  }
  else{
    return zneg_fd->GetBinContent(zneg_fd->FindBin(xval,yval,zval));
  }
}
float zoffsetbd(float xval,float yval,float zval){
  if(xval>=0){
    return zpos_bd->GetBinContent(zpos_bd->FindBin(xval,yval,zval));
  }
  else{
    return zneg_bd->GetBinContent(zneg_bd->FindBin(xval,yval,zval));
  }
}

void protoDUNE_validate_calib::Loop(int mn)
{

  if (fChain == 0) return;
  fChain->GetEntry(0);
  if (run>10000) run = 0;

  double ref_dQdx[3] = {65.75, 63.5, 59.29};
  double calib_const[3]; // = {1.162e-3, 1.116e-3, 1.0405e-3};
  //double calib_const[3] = {1.116e-3, 1.162e-3, 1.0405e-3}; //wrong
  //  if (run==0){//mc
  //    calib_const[0] = 1.036e-3;
  //    calib_const[1] = 1.034e-3;
  //    calib_const[2] = 1.0205e-3;
  //  }
  double median_dQdx[3];

  ////////////////////// Importing Y-Z plane fractional corrections /////////////

  TH2F *YZ_negativeX_corr[3];
  TH2F *YZ_positiveX_corr[3];
  TH1F *X_corr[3];
  TFile *my_file, *my_file2;

  if (recalib){
    ifstream in;
    for (int i = 0; i<3; ++i){
      in.open(Form("global_median_%d_r%d.txt",i, run));
      char line[1024];
      while(1){
        in.getline(line,1024);
        if(!in.good()) break;
        int n;
        sscanf(line, "%d %lf", &n, &median_dQdx[i]);
      }
      in.close();
      in.clear();
      in.open(Form("calconst_%d.txt", i));
      while(1){
        in.getline(line,1024);
        if(!in.good()) break;
        double error;
        sscanf(line, "%lf %lf", &calib_const[i], &error);
      }
      in.close();
      in.clear();
      cout<<"median_dQdx["<<i<<"]="<<median_dQdx[i]<<" calorimetry constant = "<<calib_const[i]<<endl;
    }

    my_file = TFile::Open(Form("YZcalo_mich%d_r%d.root",mn, run));
    my_file2 = TFile::Open(Form("Xcalo_mich%d_r%d.root",mn, run));
  
    for (int i = 0; i<3; ++i){
      YZ_negativeX_corr[i] = (TH2F*)my_file->Get(Form("correction_dqdx_ZvsY_negativeX_hist_%d",i));
      YZ_positiveX_corr[i] = (TH2F*)my_file->Get(Form("correction_dqdx_ZvsY_positiveX_hist_%d",i));
      X_corr[i] = (TH1F*)my_file2->Get(Form("dqdx_X_correction_hist_%d",i));
    }
  }
  
  //int x_bin_size=5;
  //int y_bin_size = 5; // nbiny bins in y direction
  //int z_bin_size = 5; // nbinz bins in z direction
  std::cout<<"efield at the anode neg: "<<tot_Ef(-352,300,300)<<std::endl;
  std::cout<<"efield at the anode pos: "<<tot_Ef(352,300,300)<<std::endl;

  //TFile *file = new TFile(Form("Validate_mich%d_r%d.root",mn, run),"recreate");
  TFile *file = new TFile(outfile.c_str(),"recreate");

  TTree *tree = new TTree("calotree","calo tree");
  short plane;
  float dQdx;
  float dQdx_e;
  float dEdx;
  float KE;
  float RR;
  float E;
  tree->Branch("plane",&plane,"plane/S");
  tree->Branch("dQdx",&dQdx,"dQdx/F");
  tree->Branch("dQdx_e",&dQdx_e,"dQdx_e/F");
  tree->Branch("dEdx",&dEdx,"dEdx/F");
  tree->Branch("KE",&KE,"KE/F");
  tree->Branch("RR",&RR,"RR/F");
  tree->Branch("E",&E,"E/F");

  TH1D *deltax = new TH1D("deltax","Reco_x - true_x (cm)", 100,-10,10);
  TH1D *deltay = new TH1D("deltay","Reco_y - true_y (cm)", 100,-10,10);
  TH1D *deltaz = new TH1D("deltaz","Reco_z - true_z (cm)", 100,-10,10);
  TH1D *dist = new TH1D("dist","dist",100,0,10);
  deltax->Sumw2();
  deltay->Sumw2();
  deltaz->Sumw2();
  dist->Sumw2();
  TH1D *dist_pl[3];
  for (int i = 0; i<3; ++i){
    dist_pl[i] = new TH1D(Form("dist_pl%d",i),Form("dist plane %d",i),100,-10,10);
    dist_pl[i]->Sumw2();
  }

  TH1F *dqdx_X_hist[3];
  TH1F *dedx_X_hist[3];
  TH1F *hdqdx[3];
  TH1F *hdedx[3];
  TH1F *hdqdx_mip[3];
  TH1F *hdedx_mip[3];

  vector<vector<vector<float>>> dqdx_value(3);
  vector<vector<vector<float>>> dedx_value(3);

  for (unsigned int i = 0; i<3; ++i){
    dqdx_X_hist[i] = new TH1F(Form("dqdx_X_hist_%d",i), Form("plane_%d;X Coordinate(cm);dQ/dx(ADC/cm)",i),144,-360,360);
    dedx_X_hist[i] = new TH1F(Form("dedx_X_hist_%d",i), Form("plane_%d;X Coordinate(cm);dE/dx(MeV/cm)",i),144,-360,360);


    hdqdx[i] = new TH1F(Form("hdqdx_%d",i), Form("plane_%d;dQ/dx (e/cm);Entries",i), 100,0,2e5);
    hdqdx[i]->Sumw2();
    hdedx[i] = new TH1F(Form("hdedx_%d",i), Form("plane_%d;dE/dx (MeV/cm);Entries",i), 100,0,10);
    hdedx[i]->Sumw2();
    hdqdx_mip[i] = new TH1F(Form("hdqdx_mip_%d",i), Form("plane_%d;dQ/dx (e/cm);Entries",i), 100,0,2e5);
    hdqdx_mip[i]->Sumw2();
    hdedx_mip[i] = new TH1F(Form("hdedx_mip_%d",i), Form("plane_%d;dE/dx (MeV/cm);Entries",i), 100,0,10);
    hdedx_mip[i]->Sumw2();

    dqdx_value[i].resize(144);
    dedx_value[i].resize(144);
                      
  }    

  ////////////////////////dEdx info//////////////////////////
  const size_t nbins = 50;
  double binsize = 10;

  TVectorD meanKE0(nbins);
  TVectorD meanKE1(nbins);
  TVectorD meanKE2(nbins);
  
  TVectorD meanPitch0(nbins);
  TVectorD meanPitch1(nbins);
  TVectorD meanPitch2(nbins);

  double avgKE[3][nbins] = {{0}};
  double avgPitch[3][nbins] = {{0}};
  int nhits[3][nbins] = {{0}};
  TVectorD params(2);
  params[0] = alp;
  params[1] = bet;

  std::vector<std::vector<TH1D*>> dedx(3);
  TH2D *dedxke[3];

  for (size_t i = 0; i < 3; ++i) {
    dedxke[i] = new TH2D(Form("dedxke_%zu", i),
                         Form("Plane:%zu;KE (MeV);dE/dx (MeV/cm)", i),
                         250,0,500,80,0,8);
    for (size_t j = 0; j < nbins; ++j) {
      if(j < 2) {
        dedx[i].push_back(new TH1D(Form("dedx_%zu_%zu", i, j), 
                                   Form("Plane:%zu, %d<KE<%d MeV;dE/dx (MeV/cm)", i, int(j*binsize), int((j+1)*binsize)),
                                   300, 0.0, 15));
      }
      else {
        dedx[i].push_back(new TH1D(Form("dedx_%zu_%zu", i, j),
                                   Form("Plane:%zu, %d<KE<%d MeV;dE/dx (MeV/cm)", i, int(j*binsize), int((j+1)*binsize)),
                                   200, 0.0, 10));
      }
      dedx[i].back()->Sumw2();
    }
  }
  
  const int np = 13;
  double spline_KE[np] = {10, 14, 20, 30, 40, 80, 100, 140, 200, 300, 400, 800, 1000};
  double spline_Range[np] = {0.70437, 1.27937, 2.37894, 4.72636, 7.5788, 22.0917, 30.4441, 48.2235, 76.1461, 123.567, 170.845, 353.438, 441.476};

  TSpline3 *sp = new TSpline3("Cubic Spline", spline_Range, spline_KE,13,"b2e2",0,0);

  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t real_nentries = fChain->GetEntries();
  if (real_nentries >200000){
    cout<<"Total entries = "<<real_nentries<<endl;
    cout<<"Only use 200000 events."<<endl;
    nentries = 200000;
    real_nentries = nentries;
  }
  // Long64_t nbytes = 0, nb = 0; // unused
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    // for (Long64_t jentry=0; jentry<10000;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    // nb = fChain->GetEntry(jentry);   nbytes += nb; // unused
    if(jentry%10000==0) cout<<jentry<<"/"<<real_nentries<<endl;
    
    int x_bin;
    for(int i=0; i<cross_trks; ++i){
      
      //We want to select cathode crossing tracks
      if(trkstartx[i]*trkendx[i]>=0) continue;

      bool fid_xing = (TMath::Abs(trkstartx[i])>350||trkstarty[i]<50||trkstarty[i]>550||trkstartz[i]<50||trkstartz[i]>645)&&(TMath::Abs(trkendx[i])>350||trkendy[i]<50||trkendy[i]>550||trkendz[i]<50||trkendz[i]>645);

      bool angle_2_xing = (abs(180/TMath::Pi()*trackthetaxz[i])<60 || 
                           abs(180/TMath::Pi()*trackthetaxz[i])>120) &&
        //(abs(180/TMath::Pi()*trackthetaxz[i])>10) &&
        (abs(180/TMath::Pi()*trackthetayz[i])<80 || 
         abs(180/TMath::Pi()*trackthetayz[i])>100);

      bool angle_1_neg_xing = abs(180/TMath::Pi()*trackthetaxz[i])>140;
      bool angle_1_pos_xing = abs(180/TMath::Pi()*trackthetaxz[i])<40;

      bool angle_0_neg_xing = angle_1_pos_xing;
      bool angle_0_pos_xing = angle_1_neg_xing;

      bool fid_stp = peakT_min[i]>100 && peakT_max[i]<5900 && trklen[i]>100 && trklen[i]<700 && (trkendz[i]<226 || trkendz[i]>236) && (trkstartz[i]<226 || trkstartz[i]>236) && (trkendz[i]<456 || trkendz[i]>472) && (trkstartz[i]<456 || trkstartz[i]>472);

      bool sel_stp = fid_stp && adjacent_hits[i]==0 && dist_min[i]<5;

      if (sel_stp){
        deltax->Fill(trkendx[i]-true_trkendx[i]);
        deltay->Fill(trkendy[i]-true_trkendy[i]);
        deltaz->Fill(trkendz[i]-true_trkendz[i]);
        dist->Fill(sqrt(pow(trkendx[i]-true_trkendx[i], 2)+
                        pow(trkendy[i]-true_trkendy[i], 2)+
                        pow(trkendz[i]-true_trkendz[i], 2)));
      }

      bool testneg_xing=0;
      bool testpos_xing=0;
      if(trkstartx[i]<-350 or trkendx[i]<-350) testneg_xing=1;      
      if(trkstartx[i]>350 or trkendx[i]>350) testpos_xing=1;      
      //if(!((TMath::Abs(trkstartx[i])>350||trkstarty[i]<50||trkstarty[i]>550||trkstartz[i]<50||trkstartz[i]>645)&&(TMath::Abs(trkendx[i])>350||trkendy[i]<50||trkendy[i]>550||trkendz[i]<50||trkendz[i]>645))) continue;
      
      // Loop over 3 planes
      for (short j = 0; j<3; ++j){

        plane = j;

        bool negok = true;
        bool posok = true;

        bool stpok = sel_stp;
        
        if (j == 2){
          if (!angle_2_xing){
            negok = false;
            posok = false;
          }
          stpok = stpok && lastwire[i]>5 && lastwire[i]<475;
        }
        if (j == 1){
          if (!angle_1_neg_xing) negok = false;
          if (!angle_1_pos_xing) posok = false;
          stpok = stpok && lastwire[i]>5 && lastwire[i]<795;
        }
        if (j == 0){
          if (!angle_0_neg_xing) negok = false;
          if (!angle_0_pos_xing) posok = false;
          stpok = stpok && lastwire[i]>5 && lastwire[i]<795;
        }

        double dist_corr = 0;
        if (stpok && ntrkhits[i][j]>1){
          TVector3 true_trkend(true_trkendx[i],
                               true_trkendy[i],
                               true_trkendz[i]);
          TVector3 reco_trkend(trkendx[i],
                               trkendy[i],
                               trkendz[i]);
          TVector3 trkend0, trkend1;
          if (trkresrange[i][j][0]>trkresrange[i][j][ntrkhits[i][j]-1]){
            trkend0.SetXYZ(trkhitx[i][j][ntrkhits[i][j]-1],
                           trkhity[i][j][ntrkhits[i][j]-1],
                           trkhitz[i][j][ntrkhits[i][j]-1]);
            trkend1.SetXYZ(trkhitx[i][j][ntrkhits[i][j]-2],
                             trkhity[i][j][ntrkhits[i][j]-2],
                             trkhitz[i][j][ntrkhits[i][j]-2]);
          }
          else{
            trkend0.SetXYZ(trkhitx[i][j][0],
                             trkhity[i][j][0],
                             trkhitz[i][j][0]);
            trkend1.SetXYZ(trkhitx[i][j][1],
                             trkhity[i][j][1],
                             trkhitz[i][j][1]);
          }
          if ((trkend0-true_trkend)*(trkend1-true_trkend)<0){
            dist_pl[j]->Fill((trkend0-true_trkend).Mag());
          }
          else{
            dist_pl[j]->Fill(-(trkend0-true_trkend).Mag());
          }
          if (corr_end){
            if ((trkend0-reco_trkend)*(trkend1-reco_trkend)>0){
              dist_corr = (trkend0-reco_trkend).Mag();
            }
            else{
              dist_corr = -(trkend0-reco_trkend).Mag();
            }
          }
        }
        dist_corr += rroffset[j];

        vector<float> res, dq, first5dq, last5dq;

        for(int k=0;k<ntrkhits[i][j];++k){
          trkresrange[i][j][k] += dist_corr;
          res.push_back(trkresrange[i][j][k]);
          dq.push_back(trkdqdx[i][j][k]);
        }

        if(res.size()==0) continue;

        int siz1=ntrkhits[i][j];
        float max=*max_element(res.begin(),res.end());
       
        /************************flipping wrongly ordered residual range values****************************/
        bool test=true;
        if((trkhity[i][j][siz1-1]<trkhity[i][j][0] && trkresrange[i][j][siz1-1]>trkresrange[i][j][0])||(trkhity[i][j][0]<trkhity[i][j][siz1-1] && trkresrange[i][j][0]>trkresrange[i][j][siz1-1])){
          test=false;
          for(int i1=0;i1<ntrkhits[i][j];i1++){
            trkresrange[i][j][i1]=res[siz1-i1-1];
          }
          //cout<<"This is a flipped track"<<endl;
        }

        /***************calculating the ratio of dQdx for first 5cm and last 5 cm of a track********************/ 
        for(size_t k=0;k<res.size();k++){
          if(trkresrange[i][j][k]<5) first5dq.push_back(dq[k]);
          if(trkresrange[i][j][k]>max-5) last5dq.push_back(dq[k]);
        }
       
        if(first5dq.size()<5){
          stpok = false;
        }
        else{
          float med1 = TMath::Median(first5dq.size(), &first5dq[0]);
          float med2= TMath::Median(last5dq.size(), &last5dq[0]);
          if(!((med1/med2)>1.4)) stpok = false;
        }
        if (!test) stpok = false;

        if ((fid_xing || stpok) && (negok || posok)){
          for(int k=1; k<TMath::Min(ntrkhits[i][j]-1,3000); ++k){
            if((trkhity[i][j][k]<600)&&(trkhity[i][j][k]>0)){
              if((trkhitz[i][j][k]<695)&&(trkhitz[i][j][k]>0)){
                KE = sp->Eval(trkresrange[i][j][k]);
                RR = trkresrange[i][j][k];
                if(trkhitx[i][j][k]<0 && trkhitx[i][j][k]>-360 && negok){//negative drift
                  if(trkhitx[i][j][k]<0 && trkhitx[i][j][k+1]>0) continue;
                  if(trkhitx[i][j][k]<0 && trkhitx[i][j][k-1]>0) continue;
                  float corrected_dqdx = trkdqdx[i][j][k];
                  float corrected_dedx = trkdedx[i][j][k];
                  if (recalib){
                    x_bin=X_corr[j]->FindBin(trkhitx[i][j][k]);
                  }
                  else{
                    double binwidth = (xmax-xmin)/xbins;
                    x_bin = int((trkhitx[i][j][k]-xmin)/binwidth)+1;
                    if (x_bin<1) x_bin = 1;
                    if (x_bin>xbins) x_bin = xbins;
                  }
                  if (recalib){  
                    float YZ_correction_factor_negativeX=YZ_negativeX_corr[j]->GetBinContent(YZ_negativeX_corr[j]->FindBin(trkhitz[i][j][k],trkhity[i][j][k]));
                    float X_correction_factor=X_corr[j]->GetBinContent(x_bin);
                    //float recom_correction=recom_factor(tot_Ef(trkhitx[i][j][k],trkhity[i][j][k],trkhitz[i][j][k]));
                    float normcorr = ref_dQdx[j]/median_dQdx[j];
                    if (run>10000) normcorr = 1;
                    corrected_dqdx=trkdqdx[i][j][k]*YZ_correction_factor_negativeX*X_correction_factor*normcorr/calib_const[j];
                    corrected_dedx=GetdEdx(corrected_dqdx, tot_Ef(trkhitx[i][j][k],trkhity[i][j][k],trkhitz[i][j][k]));
                  }
                  dQdx = corrected_dqdx;
                  dQdx_e = corrected_dqdx*calib_const[j];
                  dEdx = corrected_dedx;
                  E = tot_Ef(trkhitx[i][j][k],trkhity[i][j][k],trkhitz[i][j][k]);
                  if (fid_xing && testneg_xing){
                    dqdx_value[j][x_bin-1].push_back(corrected_dqdx);
                    dedx_value[j][x_bin-1].push_back(corrected_dedx);
                    hdqdx[j]->Fill(corrected_dqdx);
                    hdedx[j]->Fill(corrected_dedx);
                    //if (k==0)  cout<<event<<" neg "<<x_bin<<" "<<trkhitx[i][j][k]<<" "<<trkhity[i][j][k]<<" "<<trkhitz[i][j][k]<<endl;
                  }
                  if (stpok && trkpitch[i][j][k]>=0.5 && trkpitch[i][j][k]<=0.8){
                    size_t bin = size_t(KE/binsize);
                    //std::cout << "Bin: " << bin << std::endl;
                    if(bin < nbins){
                      dedx[j][bin]->Fill(corrected_dedx);
                      avgKE[j][bin] += KE;
                      avgPitch[j][bin] += trkpitch[i][j][k];
                      ++nhits[j][bin];
                    }
                    //dedxke[j]->Fill(trkresrange[i][j][k], corrected_dedx);
                    dedxke[j]->Fill(KE, corrected_dedx);
                    if (dEdx<100) tree->Fill();
                    if (KE>250 && KE<450){
                      hdqdx_mip[j]->Fill(corrected_dqdx);
                      hdedx_mip[j]->Fill(corrected_dedx);
                    }
                  }
                }//X containment
                if(trkhitx[i][j][k]>0 && trkhitx[i][j][k]<360 && posok){//positive drift
                  if(trkhitx[i][j][k]>0 && trkhitx[i][j][k+1]<0) continue;
                  if(trkhitx[i][j][k]>0 && trkhitx[i][j][k-1]<0) continue;
                  float corrected_dqdx = trkdqdx[i][j][k];
                  float corrected_dedx = trkdedx[i][j][k];
                  if (recalib){
                    x_bin=X_corr[j]->FindBin(trkhitx[i][j][k]);
                  }
                  else{
                    double binwidth = (xmax-xmin)/xbins;
                    x_bin = int((trkhitx[i][j][k]-xmin)/binwidth)+1;
                    if (x_bin<1) x_bin = 1;
                    if (x_bin>xbins) x_bin = xbins;
                  }
                  if (recalib){  
                    float YZ_correction_factor_positiveX=YZ_positiveX_corr[j]->GetBinContent(YZ_positiveX_corr[j]->FindBin(trkhitz[i][j][k],trkhity[i][j][k]));
                    float X_correction_factor=X_corr[j]->GetBinContent(x_bin);
                    //float recom_correction=recom_factor(tot_Ef(trkhitx[i][j][k],trkhity[i][j][k],trkhitz[i][j][k]));
                    float normcorr = ref_dQdx[j]/median_dQdx[j];
                    if (run>10000) normcorr = 1;
                    corrected_dqdx=trkdqdx[i][j][k]*YZ_correction_factor_positiveX*X_correction_factor*normcorr/calib_const[j];
                    corrected_dedx=GetdEdx(corrected_dqdx, tot_Ef(trkhitx[i][j][k],trkhity[i][j][k],trkhitz[i][j][k]));
                  }
                  dQdx = corrected_dqdx;
                  dQdx_e = corrected_dqdx*calib_const[j];
                  dEdx = corrected_dedx;
                  E = tot_Ef(trkhitx[i][j][k],trkhity[i][j][k],trkhitz[i][j][k]);
                  //if (k==2) cout<<"x = "<<trkhitx[i][j][k]<<" y = "<<trkhity[i][j][k]<<" z = "<<trkhitz[i][j][k]<<" dqdx = "<<trkdqdx[i][j][k]<<" yzcorr = "<<YZ_correction_factor_positiveX<<" xcorr = "<<X_correction_factor<<" normcorr = "<<normcorr<<" corrected_dqdx = "<<corrected_dqdx<<" corrected_dedx = "<<corrected_dedx<<" "<<YZ_positiveX_corr[j]->GetXaxis()->FindBin(trkhitz[i][j][k])<<" "<<YZ_positiveX_corr[j]->GetYaxis()->FindBin(trkhity[i][j][k])<<" "<<YZ_positiveX_corr[j]->GetBinContent(YZ_positiveX_corr[j]->GetXaxis()->FindBin(trkhitz[i][j][k]),YZ_positiveX_corr[j]->GetYaxis()->FindBin(trkhity[i][j][k]))<<" "<<x_bin<<" "<<X_correction_factor<<endl;
                  //std::cout<<TrkID[i]<<" "<<trkhitx[i][j][k]<<" "<<trkhity[i][j][k]<<" "<<trkhitz[i][j][k]<<" "<<YZ_correction_factor_positiveX_2<<" "<<trkdqdx[i][j][k]<<" "<<corrected_dqdx_2<<std::endl;
                  if (fid_xing && testpos_xing){
                    dqdx_value[j][x_bin-1].push_back(corrected_dqdx);
                    dedx_value[j][x_bin-1].push_back(corrected_dedx);
                    hdqdx[j]->Fill(corrected_dqdx);
                    hdedx[j]->Fill(corrected_dedx);
                    //if (k==0)  cout<<event<<" pos "<<x_bin<<" "<<trkhitx[i][j][k]<<" "<<trkhity[i][j][k]<<" "<<trkhitz[i][j][k]<<endl;
                  }
                  if (stpok && trkpitch[i][j][k]>=0.5 && trkpitch[i][j][k]<=0.8){
                    size_t bin = size_t(KE/binsize);
                    //std::cout << "Bin: " << bin << std::endl;
                    if(bin < nbins){
                      dedx[j][bin]->Fill(corrected_dedx);
                      avgKE[j][bin] += KE;
                      avgPitch[j][bin] += trkpitch[i][j][k];
                      ++nhits[j][bin];
                    }
                    //dedxke[j]->Fill(trkresrange[i][j][k], corrected_dedx);
                    dedxke[j]->Fill(KE, corrected_dedx);
                    if (dEdx<100) tree->Fill();
                    if (KE>250 && KE<450){
                      hdqdx_mip[j]->Fill(corrected_dqdx);
                      hdedx_mip[j]->Fill(corrected_dedx);
                    }
                  }
                }//X containment
              } // Z containment
            } // Y containment
          } // loop over hits of the track in the given plane
        }// fiducial cut
      }// loop over 3 planes
    } // loop over crossing tracks in the event
  } // loop over jentries


  ////////////////////// Getting inforamtion from dqdx_value vector ////////////////////
  for (int i = 0; i<3; ++i){
    for(size_t j=0; j<dqdx_value[i].size(); j++){
      //std::cout<<i<<" "<<j<<" "<<dqdx_value[i][k].size()<<std::endl;
      float local_median_dqdx=TMath::Median(dqdx_value[i][j].size(),&dqdx_value[i][j][0]);
      dqdx_X_hist[i]->SetBinContent(j+1,local_median_dqdx);
      float local_median_dedx=TMath::Median(dedx_value[i][j].size(),&dedx_value[i][j][0]);
      dedx_X_hist[i]->SetBinContent(j+1,local_median_dedx);
    }
    for(size_t j = 0; j<nbins; ++j){
      avgKE[i][j]/=nhits[i][j];
      avgPitch[i][j]/=nhits[i][j];
      if (i==0){
        meanKE0[j] = avgKE[i][j];
        meanPitch0[j] = avgPitch[i][j];
      }
      if (i==1){
        meanKE1[j] = avgKE[i][j];
        meanPitch1[j] = avgPitch[i][j];
      }
      if (i==2){
        meanKE2[j] = avgKE[i][j];
        meanPitch2[j] = avgPitch[i][j];
      }
    }
  }
  file->cd();
  meanKE0.Write("meanKE0");
  meanKE1.Write("meanKE1");
  meanKE2.Write("meanKE2");
  meanPitch0.Write("meanPitch0");
  meanPitch1.Write("meanPitch1");
  meanPitch2.Write("meanPitch2");
  params.Write("params");

  file->Write();
}

int main(int argc, char *argv[]) {

  bool found_fcl = false;

  std::string fcl_file;

  for (int iArg = 1; iArg < argc; iArg++) {
    if (!strcasecmp(argv[iArg],"-c")) {
     fcl_file = argv[++iArg];
     found_fcl = true;
    }
    if (!strcasecmp(argv[iArg],"-h")) {
      std::cout << "Usage: "<<argv[0] <<
                   " -c fclfile.fcl" << std::endl;
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
  outfile = pset.get<string>("outfile");
  int michelnumber = pset.get<int>("michelnumber");
  recalib = pset.get<bool>("recalib");
  corr_end = pset.get<bool>("corr_end");
  sceon = pset.get<bool>("sceon");
  xbins = pset.get<int>("xbins");
  xmin = pset.get<double>("xmin");
  xmax = pset.get<double>("xmax");
  rroffset = pset.get<vector<double>>("rroffset");
  cout<<"infile = "<<infile<<endl;
  cout<<"outfile = "<<outfile<<endl;
  cout<<"michelnumber = "<<michelnumber<<endl;
  cout<<"recalib = "<<recalib<<endl;
  cout<<"sceon = "<<sceon<<endl;
  cout<<"corr_end = "<<corr_end<<endl;
  cout<<"xbins = "<<xbins<<endl;
  cout<<"xmin = "<<xmin<<endl;
  cout<<"xmax = "<<xmax<<endl;
  for (int i = 0; i<3; ++i) cout<<"rroffset["<<i<<"] = "<<rroffset[i]<<endl;

  if (!(michelnumber == 0||michelnumber == 1||michelnumber == 2||michelnumber == 3)){
    cout << "Error: Michel tree number must be 0, 1, 2 or 3" << endl;
    return 0;
  }


  TChain* shtree = new TChain("Event");

  if (infile.substr(infile.find_last_of(".") + 1) == "root"){
    shtree->Add(Form("%s/michelremoving%d/Event", infile.c_str(), michelnumber));
  }

  else /*if(infile.substr(infile.find_last_of(".") + 1) == "list" || infile.substr(infile.find_last_of(".") + 1) == "txt")*/{
    std::ifstream in;
    in.open(infile.c_str());
    char line[1024];

    while(1){
      in.getline(line,1024);
      if (!in.good()) break;
      shtree->Add(Form("%s/michelremoving%d/Event", line, michelnumber));
    }
    in.close();
    in.clear();
  }

  protoDUNE_validate_calib t(shtree);
  
  t.Loop(michelnumber);
} // main
