#include "PDSPSystematics.h"

protoana::PDSPSystematics::PDSPSystematics(
    const std::vector<ThinSliceEvent> & events,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    std::vector<double> & beam_energy_bins,
    const std::map<std::string, ThinSliceSystematic> & pars,
    TFile & output_file) {


  SetupSyst_G4RWCoeff(pars);
  SetupSyst_BeamShift(pars, output_file);
  SetupSyst_EDivWeight(pars);
  SetupSyst_EndZNoTrackWeight(pars);
  SetupSyst_BeamMatch(pars);
  SetupSyst_BoxBeam(pars);
  SetupSyst_ELoss(pars);
  SetupSyst_ELossMuon(pars);

}

void protoana::PDSPSystematics::SetupSyst_G4RWCoeff(
    const std::map<std::string, ThinSliceSystematic> & pars) {
  for (auto it = pars.begin(); it != pars.end(); ++it) {
    if ((it->first.find("g4rw") != std::string::npos ||
         it->first.find("G4RW") != std::string::npos) &&
        (it->first.find("oeff") != std::string::npos)) {
      std::cout << "Setting up g4rw coeffsyst: " << it->first << std::endl;
    }
    else {
      continue;
    }
    fG4RWCoeffBranches[it->first] = it->second.GetOption<std::string>("Branch");
    std::cout << "Added " << it->first << " " << fG4RWCoeffBranches[it->first] << std::endl;
  }
}

double protoana::PDSPSystematics::GetSystWeight_G4RWCoeff(
    const ThinSliceEvent & event,
    const std::map<std::string, ThinSliceSystematic> & pars) {
  double weight = 1.;

  for (auto it = fG4RWCoeffBranches.begin(); it != fG4RWCoeffBranches.end();
       ++it) {
    weight *= event.GetG4RWCoeffWeight(
        it->second, pars.at(it->first).GetValue());
    //std::cout << event.GetEventID() << " " << event.GetSubrunID() << " " <<
    //             event.GetRunID() << " " <<
    //             pars.at(it->first).GetValue() << " " << weight << std::endl;
  }
  //return weight;
  return CheckAndReturn(weight, "G4RWCoeff");
}

void protoana::PDSPSystematics::SetupSyst_BeamShift(
    const std::map<std::string, ThinSliceSystematic> & pars,
    TFile & output_file) {
  if (pars.find("beam_shift") == pars.end()) {
    return;
  }

  TFile shift_file(
      pars.at("beam_shift").GetOption<std::string>("ShiftFile").c_str());
  //fSystBeamShiftMap = (TGraph2D*)shift_file.Get("g2d");
  //fSystBeamShiftMap->SetDirectory(0);
  fSystBeamShiftLimits
      = pars.at("beam_shift").GetOption<std::pair<double, double>>("Limits");

  //fSystBeamShiftRatioLimitUp = pars.at("beam_shift").GetOption<double>("RatioLimitUp");
  //fSystBeamShiftRatioLimitDown = pars.at("beam_shift").GetOption<double>("RatioLimitDown");
  fSystBeamShiftMeans = (TGraph*)shift_file.Get("gMeans");
  //fSystBeamShiftMeans->SetDirectory(0);
  fSystBeamShiftWidths = (TGraph*)shift_file.Get("gWidths");
  //fSystBeamShiftWidths->SetDirectory(0);
  shift_file.Close();

  fSystBeamShiftWeightCap = pars.at("beam_shift").GetOption<double>("WeightCap");
  /*
  fSystBeamShiftTreeSave = pars.at("beam_shift").GetOption<bool>("SaveInfo");
  if (fSystBeamShiftTreeSave) {
    output_file.mkdir("SystBeamShift");
    output_file.cd("SystBeamShift");
    fSystBeamShiftTree = new TTree("tree", "");
    fSystBeamShiftTree->Branch("Weight", &fSystBeamShiftWeight);
    fSystBeamShiftTree->Branch("Val", &fSystBeamShiftVal);
    fSystBeamShiftTree->Branch("R", &fSystBeamShiftR);
  }
  */
  fSetupSystBeamShift = true;

}

double protoana::PDSPSystematics::GetSystWeight_BeamShift(
    const ThinSliceEvent & event,
    const std::map<std::string, ThinSliceSystematic> & pars) {
  if (pars.find("beam_shift") == pars.end()) return 1.;
  if (event.GetPDG() != 211) return 1.;
  double x_val = pars.at("beam_shift").GetValue();
  double y_val = (event.GetBeamInstP() - event.GetTrueStartP())/
                  event.GetTrueStartP();
  if (y_val < fSystBeamShiftLimits.first ||
      y_val > fSystBeamShiftLimits.second) {
    return 1.;
  }

  double nominal_mean = fSystBeamShiftMeans->Eval(0.);
  double nominal_width = fSystBeamShiftWidths->Eval(0.);
  double varied_mean = fSystBeamShiftMeans->Eval(x_val);
  double varied_width = fSystBeamShiftWidths->Eval(x_val);

  double weight = (nominal_width/varied_width)*
                  exp(.5*std::pow(((y_val - nominal_mean)/nominal_width), 2)
                      - .5*std::pow(((y_val - varied_mean)/varied_width), 2));

  if (weight > fSystBeamShiftWeightCap && fSystBeamShiftWeightCap > 0.) {
    weight = fSystBeamShiftWeightCap;
  }

  /*
  if (fSystBeamShiftTreeSave) {
    fSystBeamShiftWeight = weight;
    fSystBeamShiftVal = x_val;
    fSystBeamShiftR = y_val;
    fSystBeamShiftTree->Fill();
  }*/
  //return weight;
  return CheckAndReturn(weight, "BeamShift");
}

void protoana::PDSPSystematics::SetupSyst_EDivWeight(
    const std::map<std::string, ThinSliceSystematic> & pars) {
  if (pars.find("ediv_weight") == pars.end()) {
    return;
  }
  fEDivF = pars.at("ediv_weight").GetOption<double>("F");
  fEDivCut = pars.at("ediv_weight").GetOption<double>("Cut");
}

double protoana::PDSPSystematics::GetSystWeight_EDiv(
    const ThinSliceEvent & event,
    const std::map<std::string, ThinSliceSystematic> & pars,
    int ediv_selection_ID) {
  
  if (pars.find("ediv_weight") == pars.end()) return 1.;

  const int selection_ID = event.GetSelectionID();
  if (selection_ID != ediv_selection_ID/*4*/) return 1.;

  const double endZ = event.GetRecoEndZ();
 
  double weight = 1.;
  double var = pars.at("ediv_weight").GetValue();
  if (endZ < fEDivCut) {
    weight = var;
  }
  else {
    weight = (1. - var*fEDivF)/(1. - fEDivF);
  }

  //return weight;
  return CheckAndReturn(weight, "EDivWeight");
}

void protoana::PDSPSystematics::SetupSyst_EndZNoTrackWeight(
    const std::map<std::string, ThinSliceSystematic> & pars) {
  if (pars.find("end_z_no_track_weight") != pars.end()) {
    fEndZNoTrackCut = pars.at("end_z_no_track_weight").GetOption<double>("Cut");

    std::vector<std::pair<int, std::vector<double>>> temp
        = pars.at("end_z_no_track_weight").GetOption
            <std::vector<std::pair<int, std::vector<double>>>>("Fractions");
    fEndZFractions = std::map<int, std::vector<double>>(temp.begin(), temp.end());
    std::cout << "EndZ Fracs" << std::endl;
    for (auto it = fEndZFractions.begin(); it != fEndZFractions.end(); ++it) {
      std::cout << it->first << " ";
      for (const auto & f : it->second) {
        std::cout << f << " "; 
      }
      std::cout << std::endl;
    }

  }
}

double protoana::PDSPSystematics::GetSystWeight_EndZNoTrack(
    const ThinSliceEvent & event,
    int signal_index,
    const std::map<std::string, ThinSliceSystematic> & pars,
    int upstream_ID, int no_track_ID) {
  if (pars.find("end_z_no_track_weight") == pars.end())
    return 1.;

  //Upstream Interactions
  if (event.GetSampleID() == upstream_ID) return 1.;

  if (fEndZFractions[event.GetSampleID()].size() == 0) return 1.;

  if (event.GetTrueEndZ() > fEndZNoTrackCut) return 1.;

  double variation = pars.at("end_z_no_track_weight").GetValue();

  double fraction = (signal_index > -1 ?
                     fEndZFractions[event.GetSampleID()][signal_index] :
                     fEndZFractions[event.GetSampleID()].back());

  if (fraction < 0.) {
    //std::cout << "neg fraction: " << fraction << std::endl;
    return 1.;
  }

  double weight = 
  /*return*/ (event.GetSelectionID() == no_track_ID ? variation :
          (1. - variation*fraction)/(1. - fraction));
  //for (auto it = fEndZFractions.begin(); it != fEndZFractions.end(); ++it) {
  //  std::cout << it->first << " ";
  //  for (const auto & f : it->second) {
  //    std::cout << f << " "; 
  //  }
  //  std::cout << std::endl;
  //}
  if (weight < 0.) {
    std::cout << event.GetSelectionID() << " " << event.GetSampleID() << " " <<
                 signal_index << " " << variation << " " <<
                 fraction << std::endl;
  }
  return CheckAndReturn(weight, "EndZNoTrack");

  //return pars.at("end_z_no_track_weight").GetValue();
}

void protoana::PDSPSystematics::SetupSyst_BeamMatch(
    const std::map<std::string, ThinSliceSystematic> & pars) {
  if (pars.find("beam_match_weight") == pars.end()) return; 

  fBeamMatchLimits = pars.at("beam_match_weight")
      .GetOption<std::vector<double>>("Limits");
  fBeamMatchFractions = pars.at("beam_match_weight")
      .GetOption<std::vector<double>>("Fractions");
  if (fBeamMatchFractions.size() != fBeamMatchLimits.size() + 1) {
    std::cout << "Error in Beam Match Syst Setup: Limits and Fractions don't line up" << std::endl;
    std::exception e;
    throw e;
  }
}

double protoana::PDSPSystematics::GetSystWeight_BeamMatch(
    const ThinSliceEvent & event,
    const std::map<std::string, ThinSliceSystematic> & pars,
    int upstream_ID, int no_track_ID) {
  if (pars.find("beam_match_weight") == pars.end())
    return 1.;

  //Upstream Interactions
  if (event.GetSampleID() == upstream_ID) return 1.;

  //No reco track
  if (event.GetSelectionID() == no_track_ID) return 1.;

  bool matched = (event.GetTrueID() == event.GetRecoToTrueID());
  double variation = pars.at("beam_match_weight").GetValue();

  double endz = event.GetTrueEndZ();

  int bin = -1;
  if (endz < fBeamMatchLimits[0]) {
    bin = 0;
  }
  else if (endz > fBeamMatchLimits.back()) {
    bin = fBeamMatchFractions.size() - 1;
  }
  else {
    for (size_t i = 1; i < fBeamMatchLimits.size(); ++i) {
      if (fBeamMatchLimits[i-1] < endz && endz < fBeamMatchLimits[i]) {
        bin = i;
        break;
      }
    }
  }

  double fraction = fBeamMatchFractions[bin];
  //return (matched ? variation : (1. - variation*fraction)/(1. - fraction));
  double weight = (matched ?
                   variation : (1. - variation*fraction)/(1. - fraction));
  //return (matched ? variation : (1. - variation*fraction)/(1. - fraction));
  return CheckAndReturn(weight, "BeamMatch");
}

double protoana::PDSPSystematics::GetSystWeight_UpstreamInt(
    const ThinSliceEvent & event,
    const std::map<std::string, ThinSliceSystematic> & pars,
    int upstream_ID) {
  if (pars.find("upstream_int_weight") == pars.end())
    return 1.;

  return ((event.GetSampleID() == upstream_ID) ?
          pars.at("upstream_int_weight").GetValue() :
          1.);
}

double protoana::PDSPSystematics::GetSystWeight_BGPions(
    const ThinSliceEvent & event,
    const std::map<std::string, ThinSliceSystematic> & pars,
    int past_FV_ID, int decay_ID) {
  if (pars.find("BG_pions_weight") == pars.end())
    return 1.;

  bool match = (event.GetSampleID() == past_FV_ID ||
                event.GetSampleID() == decay_ID);
  return (match ?
          pars.at("BG_pions_weight").GetValue() :
          1.);
}

void protoana::PDSPSystematics::SetupSyst_BoxBeam(
    const std::map<std::string, ThinSliceSystematic> & pars) {
  if (pars.find("box_beam_weight") == pars.end()) {
    return;
  }
  fBoxBeamRegions
      = pars.at("box_beam_weight").GetOption<
          std::vector<std::pair<double, double>>>("Regions");
  fBoxBeamFraction = pars.at("box_beam_weight").GetOption<double>("Fraction");
}

double protoana::PDSPSystematics::GetSystWeight_BoxBeam(
    const ThinSliceEvent & event,
    const std::map<std::string, ThinSliceSystematic> & pars,
    int beam_cut_ID) {
  if (pars.find("box_beam_weight") == pars.end()) return 1.;

  if (event.GetSelectionID() != beam_cut_ID) return 1.;

  double startY = event.GetRecoStartY();

  bool near_box_beam = false;
  for (const auto & region : fBoxBeamRegions) {
    if (region.first < startY && startY < region.second) {
      near_box_beam = true;
      break;
    }
  }

  double variation = pars.at("box_beam_weight").GetValue();
 // std::cout << "near/var: " << near_box_beam << "/" << variation << std::endl;
  //return (near_box_beam ? variation :
  //        (1. - variation*fBoxBeamFraction)/(1. - fBoxBeamFraction));
  return (near_box_beam ? variation : 1.);
}

void protoana::PDSPSystematics::SetupSyst_TrueBeamShift(
    const std::map<std::string, ThinSliceSystematic> & pars) {
  if (pars.find("true_beam_shift") == pars.end()) {return;}

  fTrueBeamBins = pars.at("true_beam_shift").GetOption<std::vector<double>>(
      "Bins");

  TFile file(
      pars.at("true_beam_shift").GetOption<std::string>("File").c_str(), "open");
  for (size_t i = 1; i < fTrueBeamBins.size(); ++i) {
    std::string spline_name = "spline_" + std::to_string(i);
    fTrueBeamSplines.push_back((TSpline3*)file.Get(spline_name.c_str()));
  }
}

double protoana::PDSPSystematics::GetSystWeight_TrueBeamShift(
    const ThinSliceEvent & event,
    const std::map<std::string, ThinSliceSystematic> & pars) {
  if (pars.find("true_beam_shift") == pars.end()) return 1.;
  if (event.GetPDG() != 211) return 1.;

  double true_startP = event.GetTrueStartP();
  int bin = 0;
  for (size_t i = 1; i < fTrueBeamBins.size(); ++i) {
    if (fTrueBeamBins[i-1] <= true_startP && true_startP <= fTrueBeamBins[i]) {
      bin = i-1;
      break;
    }
  }

  TSpline3 * spline = fTrueBeamSplines[bin];
  return spline->Eval(pars.at("true_beam_shift").GetValue());
}

double protoana::PDSPSystematics::GetSystWeight_ELoss(
    const ThinSliceEvent & event,
    const std::map<std::string, ThinSliceSystematic> & pars, int upstream_ID) {
  if (pars.find("eloss_weight") == pars.end()) return 1.;
  if (event.GetSampleID() == upstream_ID) return 1.;
  if (abs(event.GetPDG()) == 13) return 1.;

  double var = pars.at("eloss_weight").GetValue();
  double frac = fELossFractions[event.GetSampleID()];

  //return ((event.GetDeltaEToTPC() > fELossCut) ?
  //        var : (1. - var*frac)/(1. - frac));
  double weight = ((event.GetDeltaEToTPC() > fELossCut) ?
                   var : (1. - var*frac)/(1. - frac));
  return CheckAndReturn(weight, "ELoss");
}

void protoana::PDSPSystematics::SetupSyst_ELoss(
    const std::map<std::string, ThinSliceSystematic> & pars) {
  if (pars.find("eloss_weight") == pars.end()) {return;}
  fELossCut = pars.at("eloss_weight").GetOption<double>("Cut");
  auto temp = pars.at("eloss_weight")
      .GetOption<std::vector<std::pair<int, double>>>("Fractions");
  fELossFractions = std::map<int, double>(temp.begin(), temp.end());
}

double protoana::PDSPSystematics::GetSystWeight_ELossMuon(
    const ThinSliceEvent & event,
    const std::map<std::string, ThinSliceSystematic> & pars, int upstream_ID) {
  if (pars.find("eloss_muon_weight") == pars.end()) return 1.;
  if (event.GetSampleID() == upstream_ID) return 1.;
  if (abs(event.GetPDG()) != 13) return 1.;

  double var = pars.at("eloss_muon_weight").GetValue();
  double frac = fELossFractions[event.GetSampleID()];

  return ((event.GetDeltaEToTPC() > fELossCut) ?
          var : (1. - var*frac)/(1. - frac));
}

void protoana::PDSPSystematics::SetupSyst_ELossMuon(
    const std::map<std::string, ThinSliceSystematic> & pars) {
  if (pars.find("eloss_muon_weight") == pars.end()) {return;}
  fELossMuonCut = pars.at("eloss_muon_weight").GetOption<double>("Cut");
  auto temp = pars.at("eloss_muon_weight")
      .GetOption<std::vector<std::pair<int, double>>>("Fractions");
  fELossMuonFractions = std::map<int, double>(temp.begin(), temp.end());
}

double protoana::PDSPSystematics::CheckAndReturn(double weight,
                                                 std::string name) {
  if (weight < 0.) {
    std::cout << "Returning negative weight " << weight << " for syst " <<
                 name << std::endl;
  }
  return weight;
}
