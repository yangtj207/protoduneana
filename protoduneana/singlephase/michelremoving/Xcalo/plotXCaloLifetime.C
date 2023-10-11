#include "TH1.h"
#include "TGraph.h"
#include "TH2.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGraphErrors.h"
#include "TVectorD.h"
#include "TTimeStamp.h"
#include <fstream>
#include "TMinuit.h"
#include "TString.h"
#include <vector>
#include <string.h>
#include "TLatex.h"
#include "TPaveStats.h"
#include "TDatime.h"
#include "TColor.h"
#include "TProfile.h"
#include "TProfile2D.h"
void plotXCaloLifetime(std::string runNumber, std::string txtLabel)
{

	TCanvas c1=TCanvas();
	ofstream outfile;
	outfile.open(Form("mich1_lifetimeresults_%s.txt",txtLabel.c_str()),std::ios_base::app);
	gROOT->LoadMacro("protoDUNEStyle.C");
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetTitleX(0.35);
	gStyle->SetOptFit(111);
	gStyle->SetStatX(0.85);
	gStyle->SetStatY(0.88); 
	std::string fileName="Xcalo_mich1_r13857466.root";
	std::string set=runNumber;
	std::string s=runNumber;
	if (runNumber.find(Form("MC"))!=std::string::npos){
	std::string set="MC";
	std::string s="MC";
	fileName="Xcalo_mich1_r13857466.root";
	}
        else{  
	fileName=Form("Xcalo_mich1_r%s.root",runNumber.c_str());
	}
	TFile my_file1(Form("%s",fileName.c_str()));
	TTree* t1=(TTree*)my_file1.Get("t1");
	Int_t hour_min_sec; Int_t year_month_day;
	t1->SetBranchAddress("hour_minute_second1",&hour_min_sec);
	t1->SetBranchAddress("year_month_date1",&year_month_day);
	t1->GetEntry(0);
	std::cout<<hour_min_sec<<','<<year_month_day<<std::endl;

	TH1F *x_correction_factor=(TH1F*)my_file1.Get("dqdx_mpv_X_hist_2");
        int nBins=x_correction_factor->GetNbinsX();
	TProfile* rebinHist=new TProfile("rebinHist","rebinHist",nBins/2,0.0,2.3);
	TH1F* rebinHistL=new TH1F("rebinHistL","rebinHistL",nBins/2,0.0,2.3);
	TH1F* rebinHistR=new TH1F("rebinHistR","rebinHistR",nBins/2,0.0,2.3);
	double velocity=156;
    std::vector<double> time_BL;
    std::vector<double> time_BR;
    std::vector<double> beamLeft;
    std::vector<double> beamRight;
    
for (int i=1; i<=x_correction_factor->GetNbinsX();i++){
    
    double x_value=(abs(x_correction_factor->GetBinCenter(1))-abs(x_correction_factor->GetBinCenter(i)))/velocity;
    double position=rebinHist->FindBin(x_value);
    

     if (x_correction_factor->GetBinContent(i)<10 || x_correction_factor->GetBinContent(i)>200) continue;
    
    
    rebinHist->Fill(x_value,x_correction_factor->GetBinContent(i));
    if (x_correction_factor->GetBinCenter(i)>0) {rebinHistL->SetBinContent(position,x_correction_factor->GetBinContent(i));
rebinHistL->SetBinError(position,x_correction_factor->GetBinError(i));
    time_BL.push_back(rebinHistL->GetBinCenter(position));
    beamLeft.push_back(rebinHistL->GetBinContent(position));    
    
    }
    if (x_correction_factor->GetBinCenter(i)<0){ rebinHistR->SetBinContent(position,x_correction_factor->GetBinContent(i));
 rebinHistR->SetBinError(position,x_correction_factor->GetBinError(i));
   time_BR.push_back(rebinHistR->GetBinCenter(position));
   beamRight.push_back(rebinHistR->GetBinContent(position));
    }
    }

	rebinHistR->GetXaxis()->SetTitle("Hit Time [ms]");
	rebinHistR->GetYaxis()->SetTitle("dQ/dx [(ADC count)#timestick/cm]");
	rebinHistR->GetYaxis()->SetRangeUser(0,100);
	rebinHistR->GetXaxis()->SetRangeUser(0,2.00);
	rebinHistL->GetXaxis()->SetRangeUser(0,2.0);
	rebinHist->GetXaxis()->SetRangeUser(0,2.0);
	rebinHistR->GetXaxis()->CenterTitle();
	rebinHistR->GetYaxis()->CenterTitle(); 
	rebinHistR->SetTitle(Form("Run %s BR",s.c_str()));
	if (runNumber.find(Form("MC"))!=std::string::npos)  rebinHist->SetTitle(Form("Prod. 4 %s BR",s.c_str()));
	double lifetimeEst_BR=(time_BR.at(time_BR.size()-1)-time_BR.at(0))/(TMath::Log(beamRight.at(0)/beamRight.at(beamRight.size()-1)));
        
        double initial_BR=beamRight.at(0);
        //std::reverse(beamLeft.begin(),beamLeft.end());
        //std::reverse(time_BL.begin(),time_BL.end());
        double initial_BL=beamLeft.at(5);
        double lifetimeEst_BL=(time_BL.at(time_BL.size()-1)-time_BL.at(0))/(TMath::Log(beamLeft.at(0)/beamLeft.at(beamLeft.size()-1)));
        if(lifetimeEst_BL<0) lifetimeEst_BL=100;
        if(lifetimeEst_BR<0) lifetimeEst_BR=100;
        std::cout<<initial_BL<<","<<lifetimeEst_BL<<std::endl;
        rebinHistR->Draw("p e0");
	TF1 *fit= new TF1("exp","[0]*exp(x/(-[1]))",0.0,2.0);
	fit->SetParName(0,"Constant");
	fit->SetParName(1,"e^{-} Lifetime [ms]");
	fit->SetParameter(0,initial_BR);
	fit->SetParameter(1,lifetimeEst_BR);
	TLatex tL;
	tL.SetNDC();  
	tL.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
	tL.Draw("SAME");
	rebinHistR->Fit("exp","WQ");
	double fullLife_R=fit->GetParameter(1); double fullConst_R=fit->GetParameter(0);
	double fullSystTau_R=fit->GetParError(1); double fullConstErr_R=fit->GetParError(0);
	c1.Print(Form("fitXCalo%s_BR.png",s.c_str()));
	c1.Print(Form("fitXCalo%s_BR.pdf",s.c_str()));
	rebinHistL->GetXaxis()->SetTitle("Hit Time [ms]");
	rebinHistL->GetYaxis()->SetTitle("dQ/dx [(ADC count)#timestick/cm]");
	rebinHistL->GetYaxis()->SetRangeUser(0,100);
	rebinHistL->GetXaxis()->CenterTitle();
	rebinHistL->GetYaxis()->CenterTitle(); 
	rebinHistL->SetTitle(Form("Run %s BL",s.c_str()));
	if (runNumber.find(Form("MC"))!=std::string::npos)  rebinHist->SetTitle(Form("Prod. 4 %s BL",s.c_str()));
	rebinHistL->Draw("p e0");
	fit->SetParameter(0,initial_BL);
	fit->SetParameter(1,lifetimeEst_BL); 
	tL.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
	tL.Draw("SAME");
	rebinHistL->Fit("exp","WQ");
	double fullLife_L=fit->GetParameter(1); double fullConst_L=fit->GetParameter(0);
	double fullSystTau_L=fit->GetParError(1); double fullConstErr_L=fit->GetParError(0);
	c1.Print(Form("fitXCalo%s_BL.png",s.c_str()));
	c1.Print(Form("fitXCalo%s_BL.pdf",s.c_str()));
	rebinHist->GetXaxis()->SetTitle("Hit Time [ms]");
	rebinHist->GetYaxis()->SetTitle("dQ/dx [(ADC count)#timestick/cm]");
	rebinHist->GetYaxis()->SetRangeUser(0,100);
	rebinHist->GetXaxis()->CenterTitle();
	rebinHist->GetYaxis()->CenterTitle(); 
	rebinHist->SetTitle(Form("Run %s",s.c_str()));
	if (runNumber.find(Form("MC"))!=std::string::npos)  rebinHist->SetTitle(Form("Prod. 4 %s",s.c_str()));
	rebinHist->Draw("p e0");
	fit->SetParameter(0,initial_BL);
	fit->SetParameter(1,lifetimeEst_BL);
	tL.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
	tL.Draw("SAME"); 
	rebinHist->Fit("exp","QW");
	double fullLife=fit->GetParameter(1); double fullConst=fit->GetParameter(0);
	double fullSystTau=fit->GetParError(1); double fullConstErr=fit->GetParError(0);
	outfile<<runNumber<<','<<year_month_day<<','<<hour_min_sec<<","<<fullLife<<","<<fullSystTau<<","<<fullConst<<","<<fullConstErr<<','<<fullLife_R<<","<<fullSystTau_R<<","<<fullConst_R<<","<<fullConstErr_R<<','<<fullLife_L<<","<<fullSystTau_L<<","<<fullConst_L<<","<<fullConstErr_L<<std::endl;

	c1.Print(Form("fitXCalo%s.png",s.c_str()));
	c1.Print(Form("fitXCalo%s.pdf",s.c_str()));
}

