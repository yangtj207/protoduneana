#include "PDSPSystematics.h"
#include <chrono>

protoana::PDSPSystematics::PDSPSystematics(
    const std::vector<ThinSliceEvent> & events,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    std::vector<double> & beam_energy_bins,
    const std::map<std::string, ThinSliceSystematic> & pars,
    TFile & output_file, int upstream_ID, int no_track_ID, int decay_ID,
    int past_FV_ID, int beam_cut_ID, int past_FV_sel_ID)
    : fUpstreamID(upstream_ID), fNoTrackID(no_track_ID), fDecayID(decay_ID),
      fPastFVID(past_FV_ID), fBeamCutID(beam_cut_ID),
      fPastFVSelectionID(past_FV_sel_ID) {

  SetupSyst_G4RWCoeff(pars);
  SetupSyst_BeamShift(pars, output_file);
  SetupSyst_EDivWeight(pars);
  SetupSyst_EndZNoTrackWeight(pars);
  SetupSyst_BeamMatch(pars);
  SetupSyst_BeamMatchHigh(pars);
  SetupSyst_BeamMatchLow(pars);
  SetupSyst_BoxBeam(pars);
  SetupSyst_ELoss(pars);
  SetupSyst_ELossMuon(pars);
  SetupSyst_QuadBeamShift(pars);

  SetupSyst_BeamShiftBins(pars);

}

double protoana::PDSPSystematics::GetEventWeight(
    const ThinSliceEvent & event,
    int signal_index,
    int selection_bin,
    const std::map<std::string, ThinSliceSystematic> & pars) {
  double weight = 1.;
  //auto begin_time = std::chrono::high_resolution_clock::now();
  for (auto it = pars.begin(); it != pars.end(); ++it) {
    if (it->second.GetIsG4RWCoeff()) {
      weight *= GetSystWeight_G4RWCoeff(event, it->second);
    }
    else if (it->second.GetIsSelVar() &&
        (event.GetSelectionID() == it->second.GetSelectionID()) &&
        (selection_bin == it->second.GetSelectionBin())) {
      weight *= it->second.GetValue();
    }
    else if (it->second.GetIsTiedG4RWCoeff()) {
      weight *= GetSystWeight_TiedG4RWCoeff(event, it->second);
    }
    else if (it->first == "ediv_weight") {
      weight *= GetSystWeight_EDiv(event, it->second);
    }
    else if (it->first == "end_z_no_track_weight") {
      weight *= GetSystWeight_EndZNoTrack(event, signal_index, it->second);
    }
    else if (it->first == "beam_match_weight") {
      weight *= GetSystWeight_BeamMatch(event, signal_index, it->second);
    }
    else if (it->first == "beam_match_low_weight") {
      weight *= GetSystWeight_BeamMatchLow(event, signal_index, it->second);
    }
    else if (it->first == "beam_match_high_weight") {
      weight *= GetSystWeight_BeamMatchHigh(event, signal_index, it->second);
    }
    else if (it->first == "eloss_weight") {
      weight *= GetSystWeight_ELoss(event, signal_index, it->second);
    }
    else if (it->first == "eloss_weight_muon") {
      weight *= GetSystWeight_ELossMuon(event, it->second);
    }
    else if (it->first == "quad_beam_shift") {
      weight *= GetSystWeight_QuadBeamShift(event, it->second);
    }
  }

  //auto end_time = std::chrono::high_resolution_clock::now();
  //auto delta =
  //    std::chrono::duration_cast<std::chrono::nanoseconds>(
  //        end_time - begin_time).count();
  //std::cout << "\tWeighting took " << delta << " ns" << std::endl;

  return weight;
}

void protoana::PDSPSystematics::SetupSyst_G4RWCoeff(
    const std::map<std::string, ThinSliceSystematic> & pars) {
  for (auto it = pars.begin(); it != pars.end(); ++it) {
    if ((it->first.find("g4rw") != std::string::npos ||
         it->first.find("G4RW") != std::string::npos) &&
        (it->first.find("oeff") != std::string::npos) &&
        (it->first.find("tied") == std::string::npos)) {
      std::cout << "Setting up g4rw coeffsyst: " << it->first << std::endl;
    }
    else {
      continue;
    }
    fG4RWCoeffBranches[it->first] = it->second.GetOption<std::string>("Branch");
    std::cout << "Added " << it->first << " " << fG4RWCoeffBranches[it->first] << std::endl;
  }
}

void protoana::PDSPSystematics::SetupSyst_BeamShiftBins(
    const std::map<std::string, ThinSliceSystematic> & pars) {
  for (auto it = pars.begin(); it != pars.end(); ++it) {
    if (it->first.find("beam_shift_bin") != std::string::npos) {
      std::cout << "Setting up beam_shift_bin syst: " << it->first << std::endl;
    }
    else {
      continue;
    }

    fBeamShiftBinMeans[it->second.GetName()] = it->second.GetOption<double>("Mean");
    fBeamShiftBinSigmas[it->second.GetName()] = it->second.GetOption<double>("Sigma");
    fBeamShiftBinRangeLows[it->second.GetName()] = it->second.GetOption<double>("RangeLow");
    fBeamShiftBinRangeHighs[it->second.GetName()] = it->second.GetOption<double>("RangeHigh");
    std::cout << "Added " << it->second.GetName() << " " << std::endl;
  }
}

//optimize this -- TODO 
double protoana::PDSPSystematics::GetSystWeight_BeamShiftBin(
    const ThinSliceEvent & event,
    const ThinSliceSystematic & par) {
  if (event.GetPDG() != 211 || event.GetTrueIncidentEnergies().size() == 0)
    return 1.;

  double beam_inst_KE = sqrt(std::pow(event.GetBeamInstP()*1.e3, 2) +
                             139.57*139.57) - 139.57;
  //std::cout << beam_inst_KE << " " << fBeamShiftBinRangeLows[par.GetName()] <<
  //             " " << fBeamShiftBinRangeHighs[par.GetName()] << " " <<
  //             par.GetName() << std::endl;

  if (beam_inst_KE < fBeamShiftBinRangeLows[par.GetName()] ||
      beam_inst_KE >= fBeamShiftBinRangeHighs[par.GetName()]) {
    //std::cout << "returning 1" << std::endl;
    return 1.;
  }

  double x = (beam_inst_KE - event.GetTrueIncidentEnergies()[0]);

  double mean = fBeamShiftBinMeans[par.GetName()];
  double stddev = fBeamShiftBinSigmas[par.GetName()];
  double weight = (exp(.5*((std::pow(x - mean, 2)) -
                   (std::pow(x - mean + par.GetValue(), 2)))/(stddev*stddev)));
  //std::cout << "Returning " << par.GetValue() << " " << weight << std::endl;
  return weight;

}

double protoana::PDSPSystematics::GetSystWeight_G4RWCoeff(
    const ThinSliceEvent & event,
    const ThinSliceSystematic & par) {
    //const std::map<std::string, ThinSliceSystematic> & pars) {
  double weight = 1.;

  //for (auto it = fG4RWCoeffBranches.begin(); it != fG4RWCoeffBranches.end();
  //     ++it) {
    weight *= event.GetG4RWCoeffWeight(
        par.GetG4RWCoeffBranch(), par.GetValue());
    //std::cout << event.GetEventID() << " " << event.GetSubrunID() << " " <<
    //             event.GetRunID() << " " << it->first << " " <<
    //             pars.at(it->first).GetValue() << " " << weight << std::endl;
  //}
  return weight;
  //return CheckAndReturn(weight, "G4RWCoeff", par);
}

double protoana::PDSPSystematics::GetSystWeight_TiedG4RWCoeff(
    const ThinSliceEvent & event,
    const ThinSliceSystematic & par) {
  double weight = 1.;

  const auto & branches = par.GetTiedG4RWBranches();
  for (const auto & br : branches) {
    //std::cout << br << " " << par.GetValue() << " " <<
    //             event.GetG4RWCoeffWeight(br, par.GetValue()) << std::endl;
    weight *= event.GetG4RWCoeffWeight(
        br, par.GetValue());
  }
  return CheckAndReturn(weight, "G4RWCoeffTied", par);
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
    const ThinSliceSystematic & par) {
    //const std::map<std::string, ThinSliceSystematic> & pars) {
  //if (pars.find("beam_shift") == pars.end()) return 1.;
  if (event.GetPDG() != 211) return 1.;
  double x_val = par/*s.at("beam_shift")*/.GetValue();
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
  return CheckAndReturn(weight, "BeamShift", par);
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
    const ThinSliceSystematic & par) {
    //const std::map<std::string, ThinSliceSystematic> & pars,
    //int ediv_selection_ID) {
  
  //if (pars.find("ediv_weight") == pars.end()) return 1.;

  const int selection_ID = event.GetSelectionID();
  if (selection_ID != fPastFVSelectionID/*ediv_selection_ID4*/) return 1.;

  const double endZ = event.GetRecoEndZ();
 
  double weight = 1.;
  double var = par/*s.at("ediv_weight")*/.GetValue();
  if (endZ < fEDivCut) {
    weight = var;
  }
  else {
    weight = (1. - var*fEDivF)/(1. - fEDivF);
  }

  if (weight < 0.) {
    std::cout << endZ << " " << fEDivCut << " " << var << " " << fEDivF <<
                 std::endl;
  }

  //return weight;
  return CheckAndReturn(weight, "EDivWeight", par);
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
    const ThinSliceSystematic & par) {

  //Upstream Interactions
  if (event.GetSampleID() == fUpstreamID) return 1.;

  if (fEndZFractions[event.GetSampleID()].size() == 0) return 1.;

  if (event.GetTrueEndZ() > fEndZNoTrackCut) return 1.;

  double variation = par.GetValue();

  double fraction = GetFractionBySample(fEndZFractions, event.GetSampleID(),
                                        signal_index);
  //double fraction = (signal_index > -1 ?
  //                   fEndZFractions[event.GetSampleID()][signal_index] :
  //                   fEndZFractions[event.GetSampleID()].back());

  if (fraction < 0.) {
    return 1.;
  }

  double weight = 
      (event.GetSelectionID() == fNoTrackID ? variation :
       (1. - variation*fraction)/(1. - fraction));
  if (weight < 0.) {
    std::cout << event.GetSelectionID() << " " << event.GetSampleID() << " " <<
                 signal_index << " " << variation << " " <<
                 fraction << std::endl;
  }
  return CheckAndReturn(weight, "EndZNoTrack", par);

}

void protoana::PDSPSystematics::SetupSyst_BeamMatch(
    const std::map<std::string, ThinSliceSystematic> & pars) {
  if (pars.find("beam_match_weight") == pars.end()) return; 

  fBeamMatchLimits = pars.at("beam_match_weight")
      .GetOption<std::vector<double>>("Limits");
  //fBeamMatchFractions = pars.at("beam_match_weight")
  //    .GetOption<std::vector<double>>("Fractions");
  auto temp = pars.at("beam_match_weight")
      .GetOption<std::vector<std::pair<int, std::vector<double>>>>("Fractions");
  fBeamMatchFractions
      = std::map<int, std::vector<double>>(temp.begin(), temp.end()); 
}

double protoana::PDSPSystematics::GetSystWeight_BeamMatch(
    const ThinSliceEvent & event,
    int signal_index,
    const ThinSliceSystematic & par) {
  ////Upstream Interactions
  //if (event.GetSampleID() == fUpstreamID) return 1.;

  //No reco track
  if (event.GetSelectionID() == fNoTrackID) return 1.;

  bool matched_to_cosmic = (event.GetRecoOrigin() == 2);
  double variation = par.GetValue();

  double fraction = GetFractionBySample(fBeamMatchFractions,
                                        event.GetSampleID(),
                                        signal_index);

  //std::cout << event.GetSampleID() << " " << signal_index << " " << fraction <<
  //             std::endl;
  double weight = (matched_to_cosmic ?
                   variation : (1. - variation*fraction)/(1. - fraction));
  return CheckAndReturn(weight, "BeamMatch", par);
}

void protoana::PDSPSystematics::SetupSyst_BeamMatchLow(
    const std::map<std::string, ThinSliceSystematic> & pars) {
  if (pars.find("beam_match_low_weight") == pars.end()) return; 

  fBeamMatchLowLimit = pars.at("beam_match_low_weight")
      .GetOption<double>("Limit");
  fBeamMatchLowFraction = pars.at("beam_match_low_weight")
      .GetOption<double>("Fraction");
  auto temp = pars.at("beam_match_low_weight")
      .GetOption<std::vector<std::pair<int, std::vector<double>>>>("Fractions");
  fBeamMatchLowFractions
      = std::map<int, std::vector<double>>(temp.begin(), temp.end()); 
  fBeamMatchLowUseSingleFrac = pars.at("beam_match_low_weight")
      .GetOption<bool>("UseSingleFrac");
}

void protoana::PDSPSystematics::SetupSyst_BeamMatchHigh(
    const std::map<std::string, ThinSliceSystematic> & pars) {
  if (pars.find("beam_match_high_weight") == pars.end()) return; 

  fBeamMatchHighLimit = pars.at("beam_match_high_weight")
      .GetOption<double>("Limit");
  fBeamMatchHighFraction = pars.at("beam_match_high_weight")
      .GetOption<double>("Fraction");
  auto temp = pars.at("beam_match_high_weight")
      .GetOption<std::vector<std::pair<int, std::vector<double>>>>("Fractions");
  fBeamMatchHighFractions
      = std::map<int, std::vector<double>>(temp.begin(), temp.end()); 
  fBeamMatchHighUseSingleFrac = pars.at("beam_match_high_weight")
      .GetOption<bool>("UseSingleFrac");
}

double protoana::PDSPSystematics::GetSystWeight_BeamMatchLow(
    const ThinSliceEvent & event,
    int signal_index,
    const ThinSliceSystematic & par) {

  //Upstream Interactions
  if (event.GetSampleID() == fUpstreamID) return 1.;

  //No reco track
  if (event.GetSelectionID() == fNoTrackID) return 1.;

  //upper limit
  if (event.GetTrueEndZ() > fBeamMatchLowLimit) return 1.;

  bool matched = (event.GetTrueID() == event.GetRecoToTrueID());
  double variation = par.GetValue();

  double fraction = (fBeamMatchLowUseSingleFrac ?
                     fBeamMatchLowFraction :
                     GetFractionBySample(
                         fBeamMatchLowFractions, event.GetSampleID(),
                         signal_index));
  //std::cout << event.GetSampleID() << " " << signal_index << " " << fraction <<
  //             std::endl;
  if (fraction < 0.) {
    return 1.; 
  }
  double weight = (matched ?
                   variation : (1. - variation*fraction)/(1. - fraction));

  /*if (fraction < 0.) {
    std::cout << "Low Fraction " << fraction << " weight " << weight <<
                 " matched? " << matched << " var " << variation << std::endl;
  }*/
  return CheckAndReturn(weight, "BeamMatchLow", par);
}

double protoana::PDSPSystematics::GetSystWeight_BeamMatchHigh(
    const ThinSliceEvent & event,
    int signal_index,
    const ThinSliceSystematic & par) {

  //Upstream Interactions
  if (event.GetSampleID() == fUpstreamID) return 1.;

  //No reco track
  if (event.GetSelectionID() == fNoTrackID) return 1.;

  //upper limit
  if (event.GetTrueEndZ() < fBeamMatchHighLimit) return 1.;

  bool matched = (event.GetTrueID() == event.GetRecoToTrueID());
  double variation = par.GetValue();

  double fraction = (fBeamMatchHighUseSingleFrac ?
                     fBeamMatchHighFraction :
                     GetFractionBySample(
                         fBeamMatchHighFractions, event.GetSampleID(),
                         signal_index));
  if (fraction < 0.) {
    return 1.; 
  }
  double weight = (matched ?
                   variation : (1. - variation*fraction)/(1. - fraction));
  /*
  if (fraction < 0.) {
    std::cout << "High Fraction " << fraction << " weight " << weight <<
                 " matched? " << matched << " var " << variation << std::endl;
  }*/
  if (weight < 0.) {
    std::cout << variation << " " << fraction << " " << event.GetSampleID() <<
                 " " << signal_index << std::endl;
  }
  return CheckAndReturn(weight, "BeamMatchHigh", par);
}



double protoana::PDSPSystematics::GetSystWeight_UpstreamInt(
    const ThinSliceEvent & event,
    const ThinSliceSystematic & par) {
    //const std::map<std::string, ThinSliceSystematic> & pars,
    //int upstream_ID) {
  //if (pars.find("upstream_int_weight") == pars.end())
  //  return 1.;

  return ((event.GetSampleID() == fUpstreamID) ?
          par/*s.at("upstream_int_weight")*/.GetValue() :
          1.);
}

double protoana::PDSPSystematics::GetSystWeight_BGPions(
    const ThinSliceEvent & event,
    const ThinSliceSystematic & par) {
    //const std::map<std::string, ThinSliceSystematic> & pars,
    //int past_FV_ID, int decay_ID) {
  //if (pars.find("BG_pions_weight") == pars.end())
  //  return 1.;

  bool match = (event.GetSampleID() == fPastFVID ||
                event.GetSampleID() == fDecayID);
  return (match ?
          par/*s.at("BG_pions_weight")*/.GetValue() :
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
    const ThinSliceSystematic & par) {
    //const std::map<std::string, ThinSliceSystematic> & pars,
    //int beam_cut_ID) {
  //if (pars.find("box_beam_weight") == pars.end()) return 1.;

  if (event.GetSelectionID() != fBeamCutID) return 1.;

  double startY = event.GetRecoStartY();

  bool near_box_beam = false;
  for (const auto & region : fBoxBeamRegions) {
    if (region.first < startY && startY < region.second) {
      near_box_beam = true;
      break;
    }
  }

  double variation = par/*s.at("box_beam_weight")*/.GetValue();
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
    const ThinSliceSystematic & par) {
    //const std::map<std::string, ThinSliceSystematic> & pars) {
  //if (pars.find("true_beam_shift") == pars.end()) return 1.;
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
  return spline->Eval(par/*s.at("true_beam_shift")*/.GetValue());
}

void protoana::PDSPSystematics::SetupSyst_QuadBeamShift(
    const std::map<std::string, ThinSliceSystematic> & pars) {
  if (pars.find("quad_beam_shift") == pars.end()) {return;}

  fQuadBeamBins = pars.at("quad_beam_shift").GetOption<std::vector<double>>(
      "Bins");

  fQuadBeamMeans = pars.at("quad_beam_shift").GetOption<std::vector<double>>(
      "Means");
  fQuadBeamSigmas = pars.at("quad_beam_shift").GetOption<std::vector<double>>(
      "Sigmas");

}

double protoana::PDSPSystematics::GetSystWeight_QuadBeamShift(
    const ThinSliceEvent & event,
    const ThinSliceSystematic & par) {
  if (event.GetPDG() != 211 || event.GetTrueIncidentEnergies().size() == 0)
    return 1.;

  double beam_inst_KE = sqrt(std::pow(event.GetBeamInstP()*1.e3, 2) +
                             139.57*139.57) - 139.57;

  int bin = -1;
  for (size_t i = 1; i < fQuadBeamBins.size(); ++i) {
    if (fQuadBeamBins[i-1] <= beam_inst_KE &&
        beam_inst_KE <= fQuadBeamBins[i]) {
      bin = i-1;
      break;
    }
  }

  //std::cout << beam_inst_KE << " " << bin << " returning 1." << std::endl;
  if (bin == -1) {
    return 1.;
  }

  double x = (beam_inst_KE - event.GetTrueIncidentEnergies()[0]);

  double mean = fQuadBeamMeans[bin];
  double stddev = fQuadBeamSigmas[bin];
  double weight = (exp(.5*((std::pow(x - mean, 2)) -
                   (std::pow(x - par.GetValue()*mean, 2)))/(stddev*stddev)));
  //std::cout << "Returning " << par.GetValue() << " " << weight << std::endl;
  return weight;
}

double protoana::PDSPSystematics::GetSystWeight_ELoss(
    const ThinSliceEvent & event,
    int signal_index,
    const ThinSliceSystematic & par) {
    //const std::map<std::string, ThinSliceSystematic> & pars, int upstream_ID) {
  //if (pars.find("eloss_weight") == pars.end()) return 1.;
  if (event.GetSampleID() == fUpstreamID) return 1.;
  if (abs(event.GetPDG()) == 13) return 1.;

  double var = par/*s.at("eloss_weight")*/.GetValue();
  //double frac = fELossFractions[event.GetSampleID()];
  double frac = (signal_index > -1 ?
                 fELossFractions[event.GetSampleID()][signal_index] :
                 fELossFractions[event.GetSampleID()].back());

  double weight = ((event.GetDeltaEToTPC() > fELossCut) ?
                   var : (1. - var*frac)/(1. - frac));
  return CheckAndReturn(weight, "ELoss", par);
}

void protoana::PDSPSystematics::SetupSyst_ELoss(
    const std::map<std::string, ThinSliceSystematic> & pars) {
  if (pars.find("eloss_weight") == pars.end()) {return;}
  fELossCut = pars.at("eloss_weight").GetOption<double>("Cut");
  /*
  auto temp = pars.at("eloss_weight")
      .GetOption<std::vector<std::pair<int, double>>>("Fractions");
  fELossFractions = std::map<int, double>(temp.begin(), temp.end());
  */

  auto temp = pars.at("eloss_weight")
      .GetOption<std::vector<std::pair<int, std::vector<double>>>>("Fractions");
  fELossFractions = std::map<int, std::vector<double>>(temp.begin(), temp.end());
  for (auto it = fELossFractions.begin(); it != fELossFractions.end(); ++it) {
    std::cout << it->first << " ";
    for (auto f : it->second) std::cout << f << " ";
    std::cout << std::endl;
  }
}

double protoana::PDSPSystematics::GetSystWeight_ELossMuon(
    const ThinSliceEvent & event,
    const ThinSliceSystematic & par) {
    //const std::map<std::string, ThinSliceSystematic> & pars, int upstream_ID) {
  //if (pars.find("eloss_muon_weight") == pars.end()) return 1.;
  if (event.GetSampleID() == fUpstreamID) return 1.;
  if (abs(event.GetPDG()) != 13) return 1.;

  double var = par/*s.at("eloss_muon_weight")*/.GetValue();
  double frac = 1.;//fELossFractions[event.GetSampleID()];

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
                                                 std::string name,
                                                 const ThinSliceSystematic & par) {
  if (weight < 0.) {
    std::cout << "Returning negative weight " << weight << " for syst " <<
                 name << std::endl;
    std::cout << "\t" << par.GetValue() << " " << par.GetName() << std::endl;
  }
  return weight;
}

double protoana::PDSPSystematics::GetFractionBySample(
    const std::map<int, std::vector<double>> & fractions, int sample_ID,
    int signal_index) {
  return (signal_index > -1 ?
          fractions.at(sample_ID)[signal_index] :
          fractions.at(sample_ID).back());
}
