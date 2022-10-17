#include "ProtoDUNEdEdXFitter.h"
#include "TString.h"
#include "TVector.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "fhiclcpp/ParameterSet.h"
#include "cetlib/filepath_maker.h"


#include "protoduneana/singlephase/michelremoving/scripts/LanGausFit.h"

ProtoDUNEdEdXFitter::ProtoDUNEdEdXFitter(
    std::string fcl_file, std::string input_file, std::string output_file)
  : fOutputFile(output_file.c_str(), "recreate"),
    fInputFile(input_file.c_str(), "open") { 

  Configure(fcl_file);

  fInputTree = (TTree*)fInputFile.Get(fTreeName.c_str());
  //fInputTree->SetDirectory(0);

  for (int i = 0; i < 40; ++i) {
    std::string name = "hdEdX_dist_" + std::to_string(i);
    fHists.push_back(
        new TH1D(name.c_str(), "", (i == 0 ? 300 : 200), 0., (i == 0 ? 15 : 10.)));
    //fHists.back()->SetDirectory(0);
    fFitOutput.push_back(0x0);
  }
  DefineFitFunction();
  MakeMinimizer();

  /*
  double init = 1.e-3;
  std::cout << fFitFunction(&init) << " " << fDOF << std::endl;
  fOutputFile.cd();
  TDirectory * out1 = fOutputFile.mkdir("Try1");
  out1->cd();
  for (const auto * h : fHists) h->Write();
  for (int i = 0; i < 40; ++i) {
    if (fFitOutput[i] != 0x0)
      fFitOutput[i]->Write(("fit_" + std::to_string(i)).c_str());
  }

  init = 1.5e-3;
  std::cout << fFitFunction(&init) << " " << fDOF << std::endl;
  fOutputFile.cd();
  TDirectory * out2 = fOutputFile.mkdir("Try2");
  out2->cd();
  for (const auto * h : fHists) h->Write();
  for (int i = 0; i < 40; ++i) {
    if (fFitOutput[i] != 0x0)
      fFitOutput[i]->Write(("fit_" + std::to_string(i)).c_str());
  }*/
}


ProtoDUNEdEdXFitter::~ProtoDUNEdEdXFitter() {
  fOutputFile.Close();
}

void ProtoDUNEdEdXFitter::RunFit() {
  std::cout << "Running Fit" << std::endl;
  bool found_minimum = false;
  try {
    found_minimum = fMinimizer->Minimize();
  }
  catch (const std::exception & e) {
    std::cerr << e.what() << std::endl;
    std::cerr << "exiting safely" << std::endl;
    found_minimum = 0;
  }

  TVector output_fit_status(1);
  output_fit_status[0] = fMinimizer->Status();
  fOutputFile.cd();
  output_fit_status.Write("fit_status");

  if (!found_minimum) {
    std::cout << "Failed to find minimum" << std::endl;
    fOutputFile.cd();
    output_fit_status.Write("fit_status");

    if (fMinimizer->Status() == 3) {
      std::cout << "Edm too large: " << fMinimizer->Edm() << std::endl;
    }
  }
  else {
    std::cout << "Running Hesse" << std::endl;
    bool hesse_good = fMinimizer->Hesse();
    std::cout << "hesse good? " << hesse_good << std::endl;
    TVector output_hesse_status(1);

    fOutputFile.cd();
    output_hesse_status[0] = (hesse_good ? 1 : 0);
    output_hesse_status.Write("hesse_status");

    std::cout << fMinimizer->VariableName(0) << " " << fMinimizer->X()[0] <<
                 " " << sqrt(fMinimizer->CovMatrix(0, 0)) <<
                 std::endl;
  }

}

void ProtoDUNEdEdXFitter::MakeMinimizer() {
  fMinimizer = std::unique_ptr<ROOT::Math::Minimizer>
    (ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));

  fMinimizer->SetMaxFunctionCalls(fMaxCalls);
  fMinimizer->SetMaxIterations(fMaxIterations);
  fMinimizer->SetTolerance(fTolerance);

  fMinimizer->SetVariable(0, "C_cal", 1.e-3, 0.01);
  fMinimizer->SetVariableLimits(0, fLowerLimit, fUpperLimit);//1.04e-3, 1.05e-3);
  fMinimizer->SetFunction(fFitFunction);

}

void ProtoDUNEdEdXFitter::Configure(std::string fcl_file) {
  fhicl::ParameterSet pset;

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
  pset = fhicl::ParameterSet::make(fcl_file, lookupPolicy);

  fDensity = pset.get<double>("Density");
  fWion = pset.get<double>("Wion");
  fBetaP = pset.get<double>("BetaP");
  fAlpha = pset.get<double>("Alpha");
  fBinSize = pset.get<double>("BinSize", 5.);
  fTreeName = pset.get<std::string>("TreeName", "tree");
  std::vector<double> kinetic_energies = pset.get<std::vector<double>>("KE");
  std::vector<double> ranges = pset.get<std::vector<double>>("Range");
  fKESpline = new TSpline3("Cubic Spline", &ranges[0], &kinetic_energies[0],
                           13, "b2e2", 0, 0);
  fTolerance = pset.get<double>("Tolerance");
  fLowerLimit = pset.get<double>("LowerLimit");
  fUpperLimit = pset.get<double>("UpperLimit");
  fMaxCalls = pset.get<int>("MaxCalls");
  fMaxIterations = pset.get<int>("MaxIterations");
  fNScanSteps = pset.get<unsigned int>("NScanSteps");
  fTargetPlane = pset.get<int>("TargetPlane", 2);
}

void ProtoDUNEdEdXFitter::DefineFitFunction() {
  fFitFunction = ROOT::Math::Functor(
      [&](double const * coeffs) {
        std::cout << "New point " << coeffs[0] << std::endl;
        for (auto * f : fFitOutput) {
          if (f != 0x0) f->Delete();
        }
        fDOF = 0;
        double chi2_result = 0.;
        for (int i = 0; i < 40; ++i) {
          double range = i*fBinSize + fBinSize/2.;
          double KE = fKESpline->Eval(range);
          if (KE > 450. || KE < 250.) continue;
          fDrawCommand.Form(
              "(exp((corrected_dq_dx/%f)*(%f*%f)/(%f*efield)) - %f)*(%f*efield/%f)",
              coeffs[0], fBetaP, fWion, fDensity, fAlpha, fDensity, fBetaP);
          int nbins = (i == 0 ? 300 : 200);
          double end = (i == 0 ? 15 : 10.);
          TString hist_cmd;
          hist_cmd.Form(">>hdEdX_dist_%d(%d, 0., %f)", i, nbins, end);
          fDrawCommand += hist_cmd;
          fCut.Form("range_bin == %d && hit_plane == %d", i, fTargetPlane);
          //std::cout << fDrawCommand << std::endl;

          if (gDirectory->Get(("hdEdX_dist_" + std::to_string(i)).c_str())) {
            gDirectory->Get(("hdEdX_dist_" + std::to_string(i)).c_str())->Delete();
          }
          fInputTree->Draw(fDrawCommand, fCut);
          fHists[i] = (TH1D*)gDirectory->Get(("hdEdX_dist_" + std::to_string(i)).c_str());

          if (fHists[i]->GetEntries() < 100) continue;

          double fr[2] = {(i == 0 ? 2.2 : 1.1), (i == 0 ? 15. : 10.)};
          double sv[4] = {0., 0., 0., 0.};
          if (fHists[i]->GetMean() < 10.) {
            sv[0] = 0.1;
            sv[1] = 1.66;
            sv[2] = fHists[i]->GetEntries()*0.05;
            sv[3]=0.05;

            if (i == 0) {
              sv[0]=0.2; sv[1]=4.7; sv[2]=20; sv[3]=.01;
            }
            else if (i == 1) {
              sv[0]=0.2; sv[1]=3.0; sv[2]=10; sv[3]=.01;
            }
            else if (i == 2) {
              sv[1]=2.5;
            }
            else if (i == 3 || i == 4){
              sv[1]=2.0;
            }
          }
          else {
            sv[0] = 0.16*fHists[i]->GetRMS();
            sv[1] = 0.9*fHists[i]->GetMean();
            sv[2] = fHists[i]->GetEntries()*100;
            sv[3] = fHists[i]->GetRMS()/5.;
          }

          double pllo[4] = {.01*sv[0], .01*sv[1], .01*sv[2], .01*sv[3]};
          double plhi[4] = {100*sv[0], 100*sv[1], 100*sv[2], 100*sv[3]};

          int status;
          double chi2;
          int ndf;
          double fp[4], fpe[4];

          //std::cout << "Running langau" << " " << fHists[i]->GetEntries() <<
          //             std::endl;
          fFitOutput[i] = langaufit(fHists[i], fr, sv, pllo, plhi, fp, fpe,
                                 &chi2, &ndf, &status);
          //std::cout << fFitOutput[i]->GetNDF() << " " <<
          //             fFitOutput[i]->GetParameter(1) << " " <<
          //             fFitOutput[i]->GetParError(1) << " " <<
          //             fFitOutput[i]->GetChisquare()/fFitOutput[i]->GetNDF() <<
          //             " " << status <<
          //             std::endl;
          if (fFitOutput[i]->GetNDF() == 0 || status < 2) continue;
          if ((fFitOutput[i]->GetParError(1) > 1000.) ||
              (fFitOutput[i]->GetChisquare()/fFitOutput[i]->GetNDF() > 10.)) continue;
          double chi2_contrib = std::pow((dpdx(KE, pitchvalue, Mmu) -
                                          fFitOutput[i]->GetParameter(1)), 2);
          chi2_contrib /= std::pow(fFitOutput[i]->GetParError(1), 2);
          chi2_result += chi2_contrib;
          fDOF++;
        }
        std::cout << coeffs[0] << " " << chi2_result << " " << fDOF << std::endl;
        return chi2_result/fDOF;
      }, 1);
}

void ProtoDUNEdEdXFitter::ParameterScans() {
  double * x = new double[fNScanSteps] {};
  double * y = new double[fNScanSteps] {};
  std::cout << "\tScanning Parameter " << fMinimizer->VariableName(0) <<
               std::endl;
  bool scanned = fMinimizer->Scan(0, fNScanSteps, x, y);
  fOutputFile.cd();
  TDirectory * out = (TDirectory *)fOutputFile.mkdir("Scans");
  out->cd();
  if (scanned) {
    TGraph gr(fNScanSteps - 2, x+1, y+1);
    gr.Write(fMinimizer->VariableName(0).c_str());
  }

  delete[] x;
  delete[] y;
}

void ProtoDUNEdEdXFitter::HandScans() {
  std::vector<double> vals, chi2s;
  std::cout << "Scanning Parameter " << fMinimizer->VariableName(0) <<
               std::endl;
  for (size_t i = 0; i < fNScanSteps; ++i) {
    vals.push_back(1.e-3 + i*1.e-6);
    chi2s.push_back(fFitFunction(&vals.back()));
  }

  fOutputFile.cd();
  TDirectory * out = (TDirectory *)fOutputFile.mkdir("Scans");
  out->cd();
  TGraph gr(vals.size() - 1, &vals[1], &chi2s[1]);
  gr.Write(fMinimizer->VariableName(0).c_str());

}
