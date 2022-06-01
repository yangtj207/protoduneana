#ifndef PROTODUNEDEDXFITTER_hh
#define PROTODUNEDEDXFITTER_hh

#include <string>
#include <vector>
#include <map>
#include "TTree.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"
#include "TSpline.h"
#include "TMath.h"
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include "TString.h"

class ProtoDUNEdEdXFitter {
 public:
  ProtoDUNEdEdXFitter(std::string fcl_file,
                      std::string input_file,
                      std::string output_file);
  ~ProtoDUNEdEdXFitter();
  void MakeMinimizer();
  void DefineFitFunction();
  void RunFit();
  void ParameterScans();
  void HandScans();
 private:

  //double fCcal = 1.e-3;

  double fDensity;
  double fWion;
  double fBetaP;
  double fAlpha;
  double fBinSize;
  TSpline3 * fKESpline /*= new TSpline3("Cubic Spline",Range,KE,13,"b2e2",0,0)*/;
  double fTolerance;
  int fTargetPlane;
  double fUpperLimit;
  double fLowerLimit;
  unsigned int fNScanSteps;
  int fMaxCalls, fMaxIterations;

  TFile fOutputFile, fInputFile;
  TTree * fInputTree;
  std::string fTreeName;
  std::vector<TH1D *> fHists;
  std::vector<TF1 *> fFitOutput;
  void Configure(std::string fcl_file);
  std::unique_ptr<ROOT::Math::Minimizer> fMinimizer;
  ROOT::Math::Functor fFitFunction;
  TString fDrawCommand, fCut;
  int fDOF;

  const double Mmu = 105.658; // Mev for Mu
  const double pitchvalue=0.65;
 
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

  double dpdx(double KE,double x,double mass){
    double g = (KE + mass)/mass;
    double b = sqrt(1. - 1./(g*g));
    double epsilon = (0.307/2.)*(18/39.948)*(x*fDensity/(b*b));
    double Me = 0.51;
    double I = 188.e-6;
    double A0 = 2.*Me*(std::pow((b*g), 2))/I;
    double A1 = epsilon/I;
    double value=(1.0/x)*epsilon*((TMath::Log(A0)) + TMath::Log(A1) + 0.2 - b*b - density(b*g));
    return value;
  };
};

#endif
