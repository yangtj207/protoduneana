#include "AbsCexDriver.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TDecompChol.h"
#include "TProfile.h"
#include "Math/ProbFunc.h"
#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"
#include <ROOT/TSeq.hxx>
#include <ROOT/TTreeProcessorMT.hxx>
#include "TTreeReader.h"
#include "TEntryList.h"

#include <thread>
#include <iomanip>
#include <iostream> 
#include <numeric>
#include "valgrind/callgrind.h"

#include "ThinSliceDriverFactory.h"

DECLARE_THINSLICEDRIVER_FACTORY_NS(protoana::AbsCexDriver, protoana, AbsCexDriver)

double protoana::AbsCexDriver::fStartXDataMean = 0.;
double protoana::AbsCexDriver::fStartXDataSigma = 0.;
double protoana::AbsCexDriver::fStartYDataMean = 0.;
double protoana::AbsCexDriver::fStartYDataSigma = 0.;
double protoana::AbsCexDriver::fStartZDataMean = 0.;
double protoana::AbsCexDriver::fStartZDataSigma = 0.;
double protoana::AbsCexDriver::fStartXMCMean = 0.;
double protoana::AbsCexDriver::fStartXMCSigma = 0.;
double protoana::AbsCexDriver::fStartYMCMean = 0.;
double protoana::AbsCexDriver::fStartYMCSigma = 0.;
double protoana::AbsCexDriver::fStartZMCMean = 0.;
double protoana::AbsCexDriver::fStartZMCSigma = 0.;

std::vector<double> fG4RWVars = {};
std::vector<std::string> fG4RWBranches = {};

protoana::AbsCexDriver::~AbsCexDriver() {}

protoana::AbsCexDriver::AbsCexDriver(
    const fhicl::ParameterSet & extra_options)
    : ThinSliceDriver(extra_options),
      fEnergyFix(extra_options.get<double>("EnergyFix")),
      fDoEnergyFix(extra_options.get<bool>("DoEnergyFix")),
      fPitch(extra_options.get<double>("WirePitch")),
      fZ0(extra_options.get<double>("Z0")),
      fMultinomial(extra_options.get<bool>("Multinomial", true)),
      fEndZCut(extra_options.get<double>("EndZCut")),
      fSliceMethod(extra_options.get<std::string>("SliceMethod")),
      fInclusive(extra_options.get<bool>("Inclusive", false)),
      fERecoSelections(extra_options.get<std::vector<int>>("ERecoSelections", {})),
      fEndZSelections(extra_options.get<std::vector<int>>("EndZSelections", {})),
      fOneBinSelections(extra_options.get<std::vector<int>>("OneBinSelections", {})),
      fBeamInstPScale(extra_options.get<double>("BeamInstPScale", 1.)),
      //fMCBeamInstPShift(extra_options.get<double>("MCBeamInstPShift", 0.)),
      fRestrictBeamInstP(extra_options.get<bool>("RestrictBeamInstP", false)),
      fDebugRestrictBeamP(extra_options.get<bool>("DebugRestrictBeamP", false)),
      fVaryDataCalibration(extra_options.get<bool>("VaryDataCalibration", false)),
      fDataCalibrationFactor(extra_options.get<double>("DataCalibrationFactor", 1.)),
      fNThreads(extra_options.get<int>("NThreads", 1)),
      fExtraHistSets(extra_options.get<std::vector<fhicl::ParameterSet>>(
          "ExtraHists", {})) {
  if (fSliceMethod == "Alt") {
    fIn = new TFile("end_slices.root", "OPEN");
    fEndSlices = (TH2D*)fIn->Get("h2D")->Clone();
    for (int i = 1; i <= 735; ++i) {
      fMeans[i] = fEndSlices->ProjectionY("", i, i)->GetMean();
    }
  }
  else if (fSliceMethod == "Traj") {
    fSliceCut = extra_options.get<int>("SliceCut");
    fTrajZStart = extra_options.get<double>("TrajZStart");
  }
  else {
    fSliceCut = std::floor((fEndZCut - (fZ0 - fPitch/2.))/fPitch);
  }


  std::cout << "EReco Selections: ";
  for (const auto & s : fERecoSelections) {
    std::cout << s << " ";
  }
  std::cout << std::endl;

  std::cout << "EndZ Selections: ";
  for (const auto & s : fEndZSelections) {
    std::cout << s << " ";
  }
  std::cout << std::endl;

  std::cout << "OneBin Selections: ";
  for (const auto & s : fOneBinSelections) {
    std::cout << s << " ";
  }
  std::cout << std::endl;

  fSkipFirstLast = extra_options.get<bool>("SkipFirstLast", false);
  fBarlowBeeston = extra_options.get<bool>("BarlowBeeston", false);
  fToSkip = extra_options.get<std::vector<int>>("ToSkip", {5, 6});
  fNWorkers = extra_options.get<size_t>("NWorkers", 1);


  //fStatVar = extra_options.get<bool>("DoMCStatVar", false);
  std::vector<fhicl::ParameterSet> cov_routines
      = extra_options.get<std::vector<fhicl::ParameterSet>>(
          "CovarianceRoutines", {});

  fFakeResolution = extra_options.get<double>("FakeResolution", -999.);
  if (fFakeResolution > 0.) fStoredEnergies = std::vector<std::vector<double>>(fNWorkers);

  for (auto & routine : cov_routines) {
    std::string name = routine.get<std::string>("Name");
    fCovarianceRoutines.push_back(name);
    if (name == "BeamShift") {
      SetupBeamShiftCovRoutine(routine);
    }
  }

  if (fFakeDataRoutine == "G4RWGrid") {
    SetupFakeDataG4RW();
    FakeDataWeight = GetFakeWeight_G4RWCoeff;
  }
  else {
    FakeDataWeight = DummyWeight; 
  }

  fThrowCalibration = extra_options.get<bool>("ThrowCalibration", false);
  fVaryCalibrationFakeData = extra_options.get<bool>("VaryCalibrationFakeData",
                                                     false);
  if (fThrowCalibration) {
    fDataCalibrationFactor = fRNG.Gaus(1., fDataCalibrationFactor);
    std::cout << "New cal factor: " << fDataCalibrationFactor << std::endl;
  }
}

void protoana::AbsCexDriver::FillMCEvents(
    TTree * tree, std::vector<ThinSliceEvent> & events,
    std::vector<ThinSliceEvent> & fake_data_events,
    int & split_val, const bool & do_split, const bool & shuffle,
    int max_entries, const bool & do_fake_data) {
  std::cout << "Filling MC Events" << std::endl;

  int sample_ID, selection_ID, event, run, subrun;
  int true_beam_PDG;
  double true_beam_interactingEnergy, reco_beam_interactingEnergy;
  double true_beam_endP, true_beam_mass, true_beam_endZ;
  double reco_beam_endZ, true_beam_startP, reco_beam_startY;
  double reco_beam_startX_SCE, reco_beam_startY_SCE, reco_beam_startZ_SCE;
  double vertex_michel_score;
  int vertex_nhits;
  double beam_inst_P;
  double leading_p_costheta, leading_piplus_costheta, leading_pi0_costheta;
  std::vector<double> * reco_beam_incidentEnergies = 0x0,
                      * true_beam_incidentEnergies = 0x0,
                      * true_beam_traj_Z = 0x0,
                      * true_beam_traj_KE = 0x0,
                      * reco_daughter_track_thetas = 0x0,
                      * reco_daughter_track_scores = 0x0;
  std::vector<int> * true_beam_slices = 0x0;
  tree->SetBranchAddress("event", &event); //good
  tree->SetBranchAddress("subrun", &subrun); //good
  tree->SetBranchAddress("run", &run); //good

  tree->SetBranchAddress("true_beam_PDG", &true_beam_PDG); //good

  if (!fInclusive) { //good
    tree->SetBranchAddress("new_interaction_topology", &sample_ID);
    tree->SetBranchAddress("selection_ID", &selection_ID);
  }
  else {
    tree->SetBranchAddress("inclusive_topology", &sample_ID);
    tree->SetBranchAddress("selection_ID_inclusive", &selection_ID);
  }
  tree->SetBranchAddress("true_beam_interactingEnergy", //good
                         &true_beam_interactingEnergy);
  tree->SetBranchAddress("true_beam_endP", &true_beam_endP); //good
  tree->SetBranchAddress("true_beam_endZ", &true_beam_endZ); //good
  tree->SetBranchAddress("true_beam_mass", &true_beam_mass); //good
  tree->SetBranchAddress("reco_beam_interactingEnergy", //good
                         &reco_beam_interactingEnergy);
  tree->SetBranchAddress("reco_beam_endZ", &reco_beam_endZ); //good
  tree->SetBranchAddress("reco_beam_calo_startX", &reco_beam_startX_SCE);
  tree->SetBranchAddress("reco_beam_calo_startY", &reco_beam_startY_SCE);
  tree->SetBranchAddress("reco_beam_calo_startZ", &reco_beam_startZ_SCE);
  tree->SetBranchAddress("reco_beam_vertex_michel_score", &vertex_michel_score);
  tree->SetBranchAddress("reco_beam_vertex_nHits", &vertex_nhits);
  tree->SetBranchAddress("reco_beam_startY", &reco_beam_startY); //good
  tree->SetBranchAddress("reco_beam_incidentEnergies", //good
                         &reco_beam_incidentEnergies);
  tree->SetBranchAddress("true_beam_incidentEnergies", //good
                         &true_beam_incidentEnergies);
  tree->SetBranchAddress("true_beam_slices", //good
                         &true_beam_slices);
  tree->SetBranchAddress("true_beam_startP", &true_beam_startP); //good
  tree->SetBranchAddress("true_beam_traj_Z", &true_beam_traj_Z); //good
  tree->SetBranchAddress("true_beam_traj_KE", &true_beam_traj_KE); //good
  std::vector<double> * calibrated_dQdX = 0x0, * beam_EField = 0x0, //good
                      * track_pitch = 0x0; //good
  tree->SetBranchAddress("reco_beam_calibrated_dQdX_SCE", &calibrated_dQdX); //good
  tree->SetBranchAddress("reco_beam_EField_SCE", &beam_EField); //good
  tree->SetBranchAddress("reco_beam_TrkPitch_SCE", &track_pitch); //good
  tree->SetBranchAddress("beam_inst_P", &beam_inst_P); //good
  tree->SetBranchAddress("reco_daughter_allTrack_Theta", &reco_daughter_track_thetas); //good
  tree->SetBranchAddress("reco_daughter_PFP_trackScore_collection", //good
                            &reco_daughter_track_scores);
  std::vector<int> * true_beam_daughter_PDG = 0x0;
  tree->SetBranchAddress("true_beam_daughter_PDG", &true_beam_daughter_PDG);//good

  int true_beam_ID, reco_beam_true_byHits_ID;
  tree->SetBranchAddress("true_beam_ID", &true_beam_ID); //good
  tree->SetBranchAddress("reco_beam_true_byHits_ID", &reco_beam_true_byHits_ID); //good

  std::vector<std::vector<double>> * g4rw_full_grid_proton_coeffs = 0x0,
                                   * g4rw_full_grid_piplus_coeffs = 0x0,
                                   * g4rw_full_fine_piplus_coeffs = 0x0,
                                   * g4rw_full_grid_abscex_coeffs = 0x0,
                                   * g4rw_primary_grid_abscex_coeffs = 0x0,
                                   * g4rw_downstream_grid_piplus_coeffs = 0x0;
  tree->SetBranchAddress("g4rw_full_grid_proton_coeffs", //good
                         &g4rw_full_grid_proton_coeffs);
  tree->SetBranchAddress("g4rw_full_grid_piplus_coeffs",
                         &g4rw_full_grid_piplus_coeffs);
  tree->SetBranchAddress("g4rw_full_fine_piplus_coeffs",
                         &g4rw_full_fine_piplus_coeffs);
  tree->SetBranchAddress("g4rw_full_grid_abscex_coeffs",
                         &g4rw_full_grid_abscex_coeffs);
  tree->SetBranchAddress("g4rw_primary_grid_abscex_coeffs",
                         &g4rw_primary_grid_abscex_coeffs);
  tree->SetBranchAddress("g4rw_downstream_grid_piplus_coeffs", //good
                         &g4rw_downstream_grid_piplus_coeffs);

  std::vector<std::vector<double>> * daughter_dQdXs = 0x0,
                                   * daughter_resRanges = 0x0,
                                   * daughter_EFields = 0x0;
  tree->SetBranchAddress(
      "reco_daughter_allTrack_calibrated_dQdX_SCE", &daughter_dQdXs);//good
  tree->SetBranchAddress(
      "reco_daughter_allTrack_resRange_SCE", &daughter_resRanges);//good
  tree->SetBranchAddress(
      "reco_daughter_allTrack_EField_SCE", &daughter_EFields);//good
  bool has_pi0_shower;
  tree->SetBranchAddress("has_shower_dist_energy", &has_pi0_shower);//good 
  std::vector<double> * true_beam_daughter_startP = 0x0;
  tree->SetBranchAddress("true_beam_daughter_startP",  &true_beam_daughter_startP);//good
  tree->SetBranchAddress("leading_p_costheta", &leading_p_costheta);//good
  tree->SetBranchAddress("leading_piplus_costheta", &leading_piplus_costheta);//good
  tree->SetBranchAddress("leading_pi0_costheta", &leading_pi0_costheta);//good
  bool is_beam_scraper;
  tree->SetBranchAddress("true_beam_is_scraper", &is_beam_scraper);//good

  std::vector<double> * reco_daughter_chi2 = 0x0,
                      * reco_daughter_truncated_dEdX = 0x0;
  std::vector<int> * reco_daughter_chi2_nhits = 0x0;
  tree->SetBranchAddress("reco_daughter_allTrack_Chi2_proton",
                         &reco_daughter_chi2);
  tree->SetBranchAddress("reco_daughter_allTrack_Chi2_ndof",
                         &reco_daughter_chi2_nhits);
  tree->SetBranchAddress("reco_daughter_allTrack_truncLibo_dEdX_pos",
                         &reco_daughter_truncated_dEdX);

  int nentries = (max_entries < 0 ? tree->GetEntries() : max_entries);
  if (max_entries > tree->GetEntries()) {
      std::string message = "Requested more entries than in MC tree";
      throw std::runtime_error(message);
  }

  split_val = nentries;

  int /*events_start = 0,*/ events_end = nentries;
  int fake_start = 0/*, fake_end = nentries*/;

  std::vector<int> event_list, fake_event_list;

  if (do_split) {
    split_val = nentries/2;
    events_end = nentries/2;
    fake_start = events_end;
    std::cout << "Note: Splitting MC in half. " <<
                 split_val << "/" << nentries <<std::endl;
    if (shuffle) {
      for (int i = 0; i < nentries; ++i) {
        event_list.push_back(i);
      }
      for (int i = 0; i < events_end; ++i) {
        int r = fRNG.Integer(event_list.size());
        fake_event_list.push_back(event_list[r]);
        event_list.erase(event_list.begin()+r);
      }

      std::sort(event_list.begin(), event_list.end());
      std::sort(fake_event_list.begin(), fake_event_list.end());

      std::cout << "Shuffled " << event_list.size() << " " <<
                   fake_event_list.size() << std::endl;
      for (int & i : event_list) {
        if (std::find(fake_event_list.begin(), fake_event_list.end(), i) !=
            fake_event_list.end()) {
          std::cout << "Incorrectly shuffled" << std::endl;
        }
      }
    }
    else {
      for (int i = 0; i < events_end; ++i) {
        event_list.push_back(i);
      }
      for (int i = fake_start; i < nentries; ++i) {
        fake_event_list.push_back(i);
      }
      std::cout << "Split " << event_list.size() << " " <<
                   fake_event_list.size() << std::endl;
    }
  }
  else {
    for (int i = 0; i < events_end; ++i) {
      event_list.push_back(i);
      fake_event_list.push_back(i);
    }
  }

  /*
  TEntryList event_list, fake_event_list;
  for (int i = events_start; i < events_end; ++i) {
    event_list.Enter(i);
  }
  for (int i = fake_start; i < fake_end; ++i) {
    fake_event_list.Enter(i);
  }

  ROOT::EnableImplicitMT(fNThreads);
  ROOT::TTreeProcessorMT tp(*tree, event_list);

  auto fill_func = [&](TTreeReader & myReader) {
    TTreeReaderValue<int> runRV(myReader, "run");
    TTreeReaderValue<int> eventRV(myReader, "event");
    TTreeReaderValue<int> sample_IDRV(myReader, (!fInclusive ? 
                                                 "new_interaction_topology" :
                                                 "inclusive_topology"));

    TTreeReaderValue<int> selection_IDRV(myReader, (!fInclusive ? 
                                                    "selection_ID" :
                                                    "selection_ID_inclusive"));
    TTreeReaderValue<int> subrunRV(myReader, "subrun");
    TTreeReaderValue<int> true_beam_PDGRV(myReader, "true_beam_PDG");

    TTreeReaderValue<double> true_beam_interactingEnergyRV(myReader, "true_beam_interactingEnergy");
    TTreeReaderValue<double> reco_beam_interactingEnergyRV(myReader, "reco_beam_interactingEnergy");
    TTreeReaderValue<double> true_beam_endPRV(myReader, "true_beam_endP");
    TTreeReaderValue<double> true_beam_massRV(myReader, "true_beam_mass");
    TTreeReaderValue<double> true_beam_endZRV(myReader, "true_beam_endZ");
    TTreeReaderValue<double> reco_beam_endZRV(myReader, "reco_beam_endZ");
    TTreeReaderValue<double> true_beam_startPRV(myReader, "true_beam_startP");
    TTreeReaderValue<double> reco_beam_startYRV(myReader, "reco_beam_startY");
    TTreeReaderValue<double> beam_inst_PRV(myReader, "beam_inst_P");
    TTreeReaderValue<double> leading_p_costhetaRV(myReader, "leading_p_costheta");
    TTreeReaderValue<double> leading_piplus_costhetaRV(myReader, "leading_piplus_costheta");
    TTreeReaderValue<double> leading_pi0_costhetaRV(myReader, "leading_pi0_costheta");

    TTreeReaderValue<std::vector<double>> reco_beam_incidentEnergiesRV(myReader, "reco_beam_incidentEnergies");
    TTreeReaderValue<std::vector<double>> true_beam_incidentEnergiesRV(myReader, "true_beam_incidentEnergies");
    TTreeReaderValue<std::vector<double>> true_beam_traj_ZRV(myReader, "true_beam_traj_Z");
    TTreeReaderValue<std::vector<double>> true_beam_traj_KERV(myReader, "true_beam_traj_KE");
    TTreeReaderValue<std::vector<double>> reco_daughter_track_thetasRV(myReader, "reco_daughter_allTrack_Theta");
    TTreeReaderValue<std::vector<double>> reco_daughter_track_scoresRV(myReader, "reco_daughter_PFP_trackScore_collection");
    TTreeReaderValue<std::vector<int>> true_beam_slicesRV(myReader, "true_beam_slices");
    TTreeReaderValue<std::vector<double>> calibrated_dQdXRV(myReader, "reco_beam_calibrated_dQdX_SCE");
    TTreeReaderValue<std::vector<double>> beam_EFieldRV(myReader, "reco_beam_EField_SCE");
    TTreeReaderValue<std::vector<double>> track_pitchRV(myReader, "reco_beam_TrkPitch_SCE");
    TTreeReaderValue<std::vector<int>> true_beam_daughter_PDGRV(myReader, "true_beam_daughter_PDG");


    TTreeReaderValue<std::vector<std::vector<double>>> g4rw_full_grid_proton_coeffsRV(myReader, "g4rw_full_grid_proton_coeffs");
    TTreeReaderValue<std::vector<std::vector<double>>> g4rw_downstream_grid_piplus_coeffsRV(myReader, "g4rw_downstream_grid_piplus_coeffs");
    TTreeReaderValue<std::vector<std::vector<double>>> daughter_dQdXsRV(myReader, "reco_daughter_allTrack_calibrated_dQdX_SCE");
    TTreeReaderValue<std::vector<std::vector<double>>> daughter_resRangesRV(myReader, "reco_daughter_allTrack_resRange_SCE");
    TTreeReaderValue<std::vector<std::vector<double>>> daughter_EFieldsRV(myReader, "reco_daughter_allTrack_EField_SCE");

    TTreeReaderValue<std::vector<double>> true_beam_daughter_startPRV(myReader, "true_beam_daughter_startP");
    TTreeReaderValue<int> true_beam_IDRV(myReader, "true_beam_ID");
    TTreeReaderValue<int> reco_beam_true_byHits_IDRV(myReader, "reco_beam_true_byHits_ID");
    TTreeReaderValue<bool> has_pi0_showerRV(myReader, "has_shower_dist_energy");
    TTreeReaderValue<bool> is_beam_scraperRV(myReader, "true_beam_is_scraper");


    std::cout << "Starting thread" << std::endl;
    while (myReader.Next()) {
      auto run = *runRV;
      auto subrun = *subrunRV;
      auto event = *eventRV;

      ThinSliceEvent thin_slice_event(event, subrun, run);
      thin_slice_event.SetSampleID(*sample_IDRV);
      thin_slice_event.SetSelectionID(*selection_IDRV);
      thin_slice_event.SetTrueInteractingEnergy(*true_beam_interactingEnergyRV);
      thin_slice_event.SetRecoInteractingEnergy(*reco_beam_interactingEnergyRV);
      thin_slice_event.SetTrueEndP(*true_beam_endPRV);
      thin_slice_event.SetTrueEndZ(*true_beam_endZRV);
      thin_slice_event.SetTrueStartP(*true_beam_startPRV);
      thin_slice_event.SetTrueMass(*true_beam_massRV);
      thin_slice_event.SetRecoEndZ(*reco_beam_endZRV);
      thin_slice_event.SetRecoStartY(*reco_beam_startYRV);

      thin_slice_event.SetTrueID(*true_beam_IDRV);
      thin_slice_event.SetRecoToTrueID(*reco_beam_true_byHits_IDRV);

      thin_slice_event.SetRecoIncidentEnergies(*reco_beam_incidentEnergiesRV);
      thin_slice_event.SetTrueIncidentEnergies(*true_beam_incidentEnergiesRV);
      thin_slice_event.SetTrueTrajZ(*true_beam_traj_ZRV);
      thin_slice_event.SetTrueTrajKE(*true_beam_traj_KERV);
      thin_slice_event.SetTrueSlices(*true_beam_slicesRV);
      thin_slice_event.SetdQdXCalibrated(*calibrated_dQdXRV);
      thin_slice_event.SetEField(*beam_EFieldRV);
      thin_slice_event.SetTrackPitch(*track_pitchRV);
      thin_slice_event.SetBeamInstP(*beam_inst_PRV);
      thin_slice_event.SetPDG(*true_beam_PDGRV);
      thin_slice_event.SetRecoDaughterTrackThetas(*reco_daughter_track_thetasRV);
      thin_slice_event.SetRecoDaughterTrackScores(*reco_daughter_track_scoresRV);
      thin_slice_event.SetHasPi0Shower(*has_pi0_showerRV);
      thin_slice_event.SetLeadingPCostheta(*leading_p_costhetaRV);
      thin_slice_event.SetLeadingPiPlusCostheta(*leading_piplus_costhetaRV);
      thin_slice_event.SetLeadingPi0Costheta(*leading_pi0_costhetaRV);
      thin_slice_event.SetIsBeamScraper(*is_beam_scraperRV);
      for (size_t j = 0; j < daughter_dQdXs->size(); ++j) {
        thin_slice_event.AddRecoDaughterTrackdQdX((*daughter_dQdXsRV)[j]);
        thin_slice_event.AddRecoDaughterTrackResRange((*daughter_resRangesRV)[j]);
        thin_slice_event.AddRecoDaughterEField((*daughter_EFieldsRV)[j]);
      }

      for (size_t j = 0; j < g4rw_downstream_grid_piplus_coeffs->size(); ++j) {
        std::string name_downstream = "g4rw_downstream_grid_piplus_coeffs_" +
                                     std::to_string(j);
        thin_slice_event.MakeG4RWCoeff(name_downstream,
                                    (*g4rw_downstream_grid_piplus_coeffsRV)[j]);
      }
      thin_slice_event.MakeG4RWCoeff("g4rw_full_grid_proton_coeffs",
                                  (*g4rw_full_grid_proton_coeffsRV)[0]);

      bool found_start = false;
      for (size_t j = 1; j < true_beam_traj_ZRV->size(); ++j) {
         
        if ((*true_beam_traj_ZRV)[j] < fTrajZStart) continue;

        double delta = (*true_beam_traj_KERV)[0] - (*true_beam_traj_KERV)[j];
        thin_slice_event.SetDeltaEToTPC(delta);
        found_start = true;
        break;
      }
      if (!found_start) {
        thin_slice_event.SetDeltaEToTPC(-999.);
      }

      std::lock_guard<std::mutex> guard(fFillMutex);
      events.push_back(thin_slice_event);
    }
  };

  tp.Process(fill_func);*/


  //for (int i = events_start; i < events_end; ++i) {
  for (int & i : event_list) {
    if (!(i % 20000)) std::cout << i << "/" << split_val << std::endl;
    tree->GetEntry(i);

    //ThinSliceEvent thin_slice_event(event, subrun, run);
    events.push_back(ThinSliceEvent(event, subrun, run));
    events.back().SetSampleID(sample_ID);
    events.back().SetSelectionID(selection_ID);
    events.back().SetTrueInteractingEnergy(true_beam_interactingEnergy);
    events.back().SetRecoInteractingEnergy(reco_beam_interactingEnergy);
    events.back().SetTrueEndP(true_beam_endP);
    events.back().SetTrueEndZ(true_beam_endZ);
    events.back().SetTrueStartP(true_beam_startP);
    events.back().SetTrueMass(true_beam_mass);
    events.back().SetRecoEndZ(reco_beam_endZ);
    events.back().SetRecoStartY(reco_beam_startY);
    events.back().SetRecoStartX_SCE(reco_beam_startX_SCE);
    events.back().SetRecoStartY_SCE(reco_beam_startY_SCE);
    events.back().SetRecoStartZ_SCE(reco_beam_startZ_SCE);

    events.back().SetVertexMichelScore(vertex_michel_score);
    events.back().SetVertexNHits(vertex_nhits);

    events.back().SetTrueID(true_beam_ID);
    events.back().SetRecoToTrueID(reco_beam_true_byHits_ID);

    events.back().SetRecoIncidentEnergies(*reco_beam_incidentEnergies);
    events.back().SetTrueIncidentEnergies(*true_beam_incidentEnergies);
    events.back().SetTrueTrajZ(*true_beam_traj_Z);
    events.back().SetTrueTrajKE(*true_beam_traj_KE);
    events.back().SetTrueSlices(*true_beam_slices);
    events.back().SetdQdXCalibrated(*calibrated_dQdX);
    events.back().SetEField(*beam_EField);
    events.back().SetTrackPitch(*track_pitch);
    events.back().SetBeamInstP(beam_inst_P);
    events.back().SetPDG(true_beam_PDG);
    events.back().SetRecoDaughterTrackThetas(*reco_daughter_track_thetas);
    events.back().SetRecoDaughterTrackScores(*reco_daughter_track_scores);
    events.back().SetHasPi0Shower(has_pi0_shower);
    events.back().SetLeadingPCostheta(leading_p_costheta);
    events.back().SetLeadingPiPlusCostheta(leading_piplus_costheta);
    events.back().SetLeadingPi0Costheta(leading_pi0_costheta);
    events.back().SetIsBeamScraper(is_beam_scraper);

    for (size_t j = 0; j < daughter_dQdXs->size(); ++j) {
      events.back().AddRecoDaughterTrackdQdX((*daughter_dQdXs)[j]);
      events.back().AddRecoDaughterTrackResRange((*daughter_resRanges)[j]);
      events.back().AddRecoDaughterEField((*daughter_EFields)[j]);
      auto nhits = (*reco_daughter_chi2_nhits)[j];
      events.back().AddOneChi2PerHit(
          (nhits > 0) ?
          (*reco_daughter_chi2)[j]/nhits :
          -999.
      );
      events.back().AddOneTrunc_dEdX(
        (*reco_daughter_truncated_dEdX)[j]
      );
    }

    for (size_t j = 0; j < g4rw_downstream_grid_piplus_coeffs->size(); ++j) {
      std::string name_downstream = "g4rw_downstream_grid_piplus_coeffs_" +
                                   std::to_string(j);
      events.back().MakeG4RWCoeff(name_downstream,
                                  (*g4rw_downstream_grid_piplus_coeffs)[j]);

      std::string name_full = "g4rw_full_grid_piplus_coeffs_" + std::to_string(j);
      events.back().MakeG4RWCoeff(name_full, (*g4rw_full_grid_piplus_coeffs)[j]);

      //std::string fine_full = "g4rw_full_fine_piplus_coeffs_" + std::to_string(j);
      //events.back().MakeG4RWCoeff(fine_full, (*g4rw_full_fine_piplus_coeffs)[j]);
    }

    for (size_t j = 0; j < g4rw_full_fine_piplus_coeffs->size(); ++j) {
      std::string fine_full = "g4rw_full_fine_piplus_coeffs_" + std::to_string(j);
      events.back().MakeG4RWCoeff(fine_full, (*g4rw_full_fine_piplus_coeffs)[j]);
    }

    for (size_t j = 0; j < g4rw_full_grid_abscex_coeffs->size(); ++j) {
      std::string abscex = "g4rw_full_grid_abscex_coeffs_" + std::to_string(j);
      events.back().MakeG4RWCoeff(abscex, (*g4rw_full_grid_abscex_coeffs)[j]);
    }

    for (size_t j = 0; j < g4rw_primary_grid_abscex_coeffs->size(); ++j) {
      std::string abscex = "g4rw_primary_grid_abscex_coeffs_" + std::to_string(j);
      events.back().MakeG4RWCoeff(abscex, (*g4rw_primary_grid_abscex_coeffs)[j]);
    }

    events.back().MakeG4RWCoeff("g4rw_full_grid_proton_coeffs",
                                (*g4rw_full_grid_proton_coeffs)[0]);

    bool found_start = false;
    for (size_t j = 1; j < true_beam_traj_Z->size(); ++j) {
       
      if ((*true_beam_traj_Z)[j] < fTrajZStart) continue;

      double delta = (*true_beam_traj_KE)[0] - (*true_beam_traj_KE)[j];
      events.back().SetDeltaEToTPC(delta);
      found_start = true;
      break;
    }
    if (!found_start) {
      events.back().SetDeltaEToTPC(-999.);
    }

    //events.push_back(thin_slice_event);
    //events.emplace_back(std::move(thin_slice_event));
  }

  if (do_fake_data) {
    std::cout << "Filling fake data" << std::endl;
    //for (int i = fake_start; i < fake_end; ++i) {
    auto start_time = std::chrono::high_resolution_clock::now();
    for (int & i : fake_event_list) {
      //std::cout << "Fake event " << i << std::endl;
      tree->GetEntry(i);

      fake_data_events.push_back(ThinSliceEvent(event, subrun, run));
      if (!(fake_data_events.size() % 20000))
        std::cout << fake_data_events.size() << std::endl;
      fake_data_events.back().SetSampleID(sample_ID);
      fake_data_events.back().SetSelectionID(selection_ID);
      fake_data_events.back().SetTrueInteractingEnergy(true_beam_interactingEnergy);
      fake_data_events.back().SetRecoInteractingEnergy(reco_beam_interactingEnergy);
      fake_data_events.back().SetTrueEndP(true_beam_endP);
      fake_data_events.back().SetTrueEndZ(true_beam_endZ);
      fake_data_events.back().SetTrueStartP(true_beam_startP);
      fake_data_events.back().SetTrueMass(true_beam_mass);
      fake_data_events.back().SetRecoEndZ(reco_beam_endZ);
      fake_data_events.back().SetRecoStartX_SCE(reco_beam_startX_SCE);
      fake_data_events.back().SetRecoStartY_SCE(reco_beam_startY_SCE);
      fake_data_events.back().SetRecoStartZ_SCE(reco_beam_startZ_SCE);
      fake_data_events.back().SetVertexMichelScore(vertex_michel_score);
      fake_data_events.back().SetVertexNHits(vertex_nhits);

      fake_data_events.back().SetRecoIncidentEnergies(*reco_beam_incidentEnergies);
      fake_data_events.back().SetTrueIncidentEnergies(*true_beam_incidentEnergies);
      fake_data_events.back().SetTrueTrajZ(*true_beam_traj_Z);
      fake_data_events.back().SetTrueTrajKE(*true_beam_traj_KE);
      fake_data_events.back().SetTrueSlices(*true_beam_slices);
      fake_data_events.back().SetdQdXCalibrated(*calibrated_dQdX);
      fake_data_events.back().SetEField(*beam_EField);
      fake_data_events.back().SetTrackPitch(*track_pitch);
      fake_data_events.back().SetBeamInstP(beam_inst_P);
      fake_data_events.back().SetPDG(true_beam_PDG);
      fake_data_events.back().SetRecoDaughterTrackThetas(*reco_daughter_track_thetas);
      fake_data_events.back().SetTrueDaughterPDGs(*true_beam_daughter_PDG);
      fake_data_events.back().SetTrueDaughterStartPs(*true_beam_daughter_startP);
      fake_data_events.back().SetRecoDaughterTrackScores(*reco_daughter_track_scores);
      fake_data_events.back().SetHasPi0Shower(has_pi0_shower);
      fake_data_events.back().SetLeadingPCostheta(leading_p_costheta);
      fake_data_events.back().SetLeadingPiPlusCostheta(leading_piplus_costheta);
      fake_data_events.back().SetLeadingPi0Costheta(leading_pi0_costheta);
      fake_data_events.back().SetIsBeamScraper(is_beam_scraper);

      fake_data_events.back().SetTrueID(true_beam_ID);
      fake_data_events.back().SetRecoToTrueID(reco_beam_true_byHits_ID);

      //fake_data_events.back().MakeG4RWBranch("g4rw_alt_primary_plus_sigma_weight",
      //                              *g4rw_alt_primary_plus_sigma_weight);
      //fake_data_events.back().MakeG4RWBranch("g4rw_alt_primary_minus_sigma_weight",
      //                              *g4rw_alt_primary_minus_sigma_weight);
      //fake_data_events.back().MakeG4RWBranch("g4rw_full_primary_plus_sigma_weight",
      //                              *g4rw_full_primary_plus_sigma_weight);
      //fake_data_events.back().MakeG4RWBranch("g4rw_full_primary_minus_sigma_weight",
      //                              *g4rw_full_primary_minus_sigma_weight);
      for (size_t j = 0; j < daughter_dQdXs->size(); ++j) {
        fake_data_events.back().AddRecoDaughterTrackdQdX((*daughter_dQdXs)[j]);
        fake_data_events.back().AddRecoDaughterTrackResRange((*daughter_resRanges)[j]);
        fake_data_events.back().AddRecoDaughterEField((*daughter_EFields)[j]);
        auto nhits = (*reco_daughter_chi2_nhits)[j];
        fake_data_events.back().AddOneChi2PerHit(
            (nhits > 0) ?
            (*reco_daughter_chi2)[j]/nhits :
            -999.
        );
        fake_data_events.back().AddOneTrunc_dEdX(
          (*reco_daughter_truncated_dEdX)[j]
        );
      }

      for (size_t j = 0; j < g4rw_full_fine_piplus_coeffs->size(); ++j) {
        std::string fine_full = "g4rw_full_fine_piplus_coeffs_" + std::to_string(j);
        fake_data_events.back().MakeG4RWCoeff(fine_full, (*g4rw_full_fine_piplus_coeffs)[j]);
      }
      for (size_t j = 0; j < g4rw_full_grid_abscex_coeffs->size(); ++j) {
        std::string abscex = "g4rw_full_grid_abscex_coeffs_" + std::to_string(j);
        fake_data_events.back().MakeG4RWCoeff(abscex, (*g4rw_full_grid_abscex_coeffs)[j]);

      }
      for (size_t j = 0; j < g4rw_primary_grid_abscex_coeffs->size(); ++j) {
        std::string abscex = "g4rw_primary_grid_abscex_coeffs_" + std::to_string(j);
        fake_data_events.back().MakeG4RWCoeff(abscex, (*g4rw_primary_grid_abscex_coeffs)[j]);
      }

      for (size_t j = 0; j < g4rw_downstream_grid_piplus_coeffs->size(); ++j) {
        //std::string name_full = "g4rw_full_grid_weights_" + std::to_string(j);
        //fake_data_events.back().MakeG4RWBranch(name_full, (*g4rw_full_grid_weights)[j]);
        std::string name_full = "g4rw_full_grid_piplus_coeffs_" + std::to_string(j);
        fake_data_events.back().MakeG4RWCoeff(name_full, (*g4rw_full_grid_piplus_coeffs)[j]);

        //std::string name_primary = "g4rw_primary_grid_weights_" +
        //                           std::to_string(j);
        //std::cout << "Adding " << name_primary << std::endl;
        //if (!(*g4rw_primary_grid_weights)[j].size())
        //  std::cout << "Adding empty branch " << event << " " << run << " " << subrun << std::endl;
        //std::string name_downstream_branch = "g4rw_downstream_grid_piplus_weights_" +
        //                             std::to_string(j);
        //fake_data_events.back().MakeG4RWBranch(name_downstream_branch,
        //                            (*g4rw_downstream_grid_piplus_weights)[j]);


        //fake_data_events.back().MakeG4RWBranch(name_primary,
        //                              (*g4rw_primary_grid_weights)[j]);
        std::string name_downstream = "g4rw_downstream_grid_piplus_coeffs_" +
                                     std::to_string(j);
        fake_data_events.back().MakeG4RWCoeff(name_downstream,
                                    (*g4rw_downstream_grid_piplus_coeffs)[j]);
      }
      //fake_data_events.back().MakeG4RWBranch("g4rw_full_grid_proton_weights",
      //                              (*g4rw_full_grid_proton_weights)[0]);
      fake_data_events.back().MakeG4RWCoeff("g4rw_full_grid_proton_coeffs",
                                            (*g4rw_full_grid_proton_coeffs)[0]);
      bool found_start = false;
      for (size_t j = 1; j < true_beam_traj_Z->size(); ++j) {
        if ((*true_beam_traj_Z)[j-1] < fTrajZStart &&
            fTrajZStart < (*true_beam_traj_Z)[j+1]) {
          double delta = (*true_beam_traj_KE)[0] - (*true_beam_traj_KE)[j];
          fake_data_events.back().SetDeltaEToTPC(delta);
          found_start = true;
          break;
        }
      }
      if (!found_start) fake_data_events.back().SetDeltaEToTPC(-999.);
    }
    auto new_time = std::chrono::high_resolution_clock::now();
    auto delta =
        std::chrono::duration_cast<std::chrono::milliseconds>(
            new_time - start_time).count();
    std::cout << "Filling fake data took " << delta << std::endl;
  }

  std::cout << "Filled MC Events" << std::endl;


}

void protoana::AbsCexDriver::BuildMCSamples(
    const std::vector<ThinSliceEvent> & events,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    std::map<int, double> & nominal_fluxes,
    std::map<int, std::vector<std::vector<double>>> & fluxes_by_sample,
    std::vector<double> & beam_energy_bins, bool use_beam_inst_P) {

  for (size_t i = 0; i < events.size(); ++i) {
    const ThinSliceEvent & event = events.at(i);

    int sample_ID = event.GetSampleID();
    int selection_ID = event.GetSelectionID();

    double true_beam_interactingEnergy = event.GetTrueInteractingEnergy();
    double reco_beam_interactingEnergy = event.GetRecoInteractingEnergy();
    double true_beam_endP = event.GetTrueEndP();
    double reco_beam_endZ = event.GetRecoEndZ();
    double true_beam_startP = event.GetTrueStartP();
    double beam_inst_P = event.GetBeamInstP();

    const std::vector<double> & reco_beam_incidentEnergies
        = event.GetRecoIncidentEnergies();
    const std::vector<double> & true_beam_incidentEnergies
        = event.GetTrueIncidentEnergies();
    const std::vector<double> & true_beam_traj_Z
        = event.GetTrueTrajZ();
    const std::vector<double> & true_beam_traj_KE
        = event.GetTrueTrajKE();
    const std::vector<int> & true_beam_slices
        = event.GetTrueSlices();

    //double beam_shift_delta = 0.;
    //if (fBeamShiftCovRoutineActive) {
    //  beam_shift_delta = GetBeamShiftDelta(true_beam_incidentEnergies);
    //  beam_inst_P += beam_shift_delta;
    //}

    double end_energy = true_beam_interactingEnergy;
    if (fSliceMethod == "Traj") {
      end_energy = sqrt(true_beam_endP*true_beam_endP*1.e6 + 139.57*139.57) - 139.57;
    }
    else if (fSliceMethod == "E") {
      end_energy = sqrt(true_beam_endP*true_beam_endP*1.e6 + 139.57*139.57) - 139.57;
    }
    else if (fSliceMethod == "Alt") {
      int bin = fEndSlices->GetXaxis()->FindBin(true_beam_traj_Z./*->*/back());
      if (bin > 0)
        end_energy = fMeans.at(bin);
    }

    if (samples.find(sample_ID) == samples.end()) {
      std::cout << "Warning: skipping sample " << sample_ID << std::endl;
      continue;
    }

    //Build the true incident energy vector based on the slice cut
    std::vector<double> good_true_incEnergies = MakeTrueIncidentEnergies(
        true_beam_traj_Z, true_beam_traj_KE, true_beam_slices,
        true_beam_incidentEnergies);
    //here
    int bin = GetBeamBin(beam_energy_bins, (use_beam_inst_P ? beam_inst_P :
                                            true_beam_startP),
                         fRestrictBeamInstP);
    if (fRestrictBeamInstP && bin == -1) continue;

    std::vector<ThinSliceSample> & samples_vec = samples.at(sample_ID)[bin];
    bool is_signal = signal_sample_checks.at(sample_ID);

    ThinSliceSample * this_sample = 0x0;
    if (!is_signal) {
      this_sample = &samples_vec.at(0);

      fluxes_by_sample[sample_ID][bin][0] += 1.;
    }
    else {
      //Iterate through the true bins and find the correct one

      bool found = false;
      if (samples_vec.size() == 1) {
        this_sample = &(samples_vec)[0];

        fluxes_by_sample[sample_ID][bin][0] += 1.;
        found = true;
      }
      else {
        for (size_t j = 1; j < samples_vec.size()-1; ++j) {
          ThinSliceSample & sample = samples_vec.at(j);
          if (sample.CheckInSignalRange(end_energy)) {
            this_sample = &sample;

            fluxes_by_sample[sample_ID][bin][j] += 1.;
            found = true;
            break;
          }
        }
      }
      if (!found) {
        //over/underflow here
        if (end_energy <= samples_vec[1].RangeLowEnd()) {
          this_sample = &samples_vec[0];

          fluxes_by_sample[sample_ID][bin][0] += 1.;
        }
        else if (end_energy >
                 samples_vec[samples_vec.size()-2].RangeHighEnd()) {
          this_sample = &samples_vec.back();
          fluxes_by_sample[sample_ID][bin].back() += 1.;
        }
        else {
          std::cout << "Warning: could not find true bin " <<
                       end_energy << std::endl;
        }
      }
    }

    /*if (sample_ID == 4) {
      std::cout << sample_ID << " " << this_sample << " " << is_signal << std::endl;
      std::cout << samples_vec.size() << std::endl;
    }*/

    int flux_type = this_sample->GetFluxType();
    nominal_fluxes[flux_type] += 1.;
    this_sample->AddFlux();
    double val[1] = {0};
    //if (selection_ID == 4) {
    if (std::find(fEndZSelections.begin(), fEndZSelections.end(), selection_ID)
        != fEndZSelections.end()) {
      TH1D * selected_hist
          = (TH1D*)this_sample->GetSelectionHists().at(selection_ID);
      if (selected_hist->FindBin(reco_beam_endZ) == 0) {
        val[0] = selected_hist->GetBinCenter(1);
      }
      else if (selected_hist->FindBin(reco_beam_endZ) >
               selected_hist->GetNbinsX()) {
        val[0] = selected_hist->GetBinCenter(selected_hist->GetNbinsX());
      }
      else {
        val[0] = reco_beam_endZ;
      }
    }
    //else if (selection_ID > 4) {
    else if (
        std::find(fOneBinSelections.begin(), fOneBinSelections.end(), selection_ID)
        != fOneBinSelections.end()) {
      val[0] = .5;
    }
    else if (reco_beam_incidentEnergies./*->*/size()) {
      double energy[1] = {reco_beam_interactingEnergy};
      if (fDoEnergyFix) {
        for (size_t k = 1; k < reco_beam_incidentEnergies.size(); ++k) {
          double deltaE = ((reco_beam_incidentEnergies)[k-1] -
                           (reco_beam_incidentEnergies)[k]);
          if (deltaE > fEnergyFix) {
            energy[0] += deltaE; 
          }
          //energy[0] += 1.e3*beam_shift_delta;
        }
      }

      TH1D * selected_hist
          = (TH1D*)this_sample->GetSelectionHists().at(selection_ID);
      if (selected_hist->FindBin(energy[0]) == 0) {
        val[0] = selected_hist->GetBinCenter(1);
      }
      else if (selected_hist->FindBin(energy[0]) >
               selected_hist->GetNbinsX()) {
        val[0] = selected_hist->GetBinCenter(selected_hist->GetNbinsX());
      }
      else {
        val[0] = energy[0];
      }
    }
    else {
      TH1D * selected_hist
          = (TH1D*)this_sample->GetSelectionHists().at(selection_ID);
      val[0] = selected_hist->GetBinCenter(1);
    }
    this_sample->FillSelectionHist(selection_ID, val);

    //Fill the total incident hist with truth info
    //this_sample->FillTrueIncidentHist(good_true_incEnergies);
    this_sample->AddIncidentEnergies(good_true_incEnergies);
    if (true_beam_incidentEnergies.size() > 0) {
    this_sample->AddESliceEnergies(
        {(true_beam_incidentEnergies)[0],
         end_energy});
    }

    FillExtraHistsMC(*this_sample, event, 1.);

    this_sample->CalculateStatVariations();

  }

}

void protoana::AbsCexDriver::RefillSampleLoop(
    const std::vector<ThinSliceEvent> & events,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    std::vector<double> & beam_energy_bins,
    const std::map<int, std::vector<double>> & signal_pars,
    const std::map<int, double> & flux_pars,
    const std::map<std::string, ThinSliceSystematic> & syst_pars,
    const std::map<std::string, ThinSliceSystematic> & g4rw_pars,
    bool fit_under_over, bool tie_under_over, bool use_beam_inst_P,
    bool fill_incident, std::map<int, TH1 *> * fix_factors,
    size_t worker_id, std::vector<size_t> n_events) {

  size_t start_event = 0;
  for (size_t i = 0; i < worker_id; ++i) {
    start_event += n_events[i];
  }

  size_t end_event = start_event + n_events[worker_id];
  size_t energy_index = 0;
  for (size_t i = start_event; i < end_event; ++i) {
    const ThinSliceEvent & event = events.at(i);
    int sample_ID = event.GetSampleID();
    int selection_ID = event.GetSelectionID();

    double true_beam_interactingEnergy = event.GetTrueInteractingEnergy();
    double reco_beam_interactingEnergy = event.GetRecoInteractingEnergy();
    double true_beam_endP = event.GetTrueEndP();
    double reco_beam_endZ = event.GetRecoEndZ();
    double true_beam_startP = event.GetTrueStartP();

    const std::vector<double> & reco_beam_incidentEnergies
        = event.GetRecoIncidentEnergies();
    const std::vector<double> & true_beam_incidentEnergies
        = event.GetTrueIncidentEnergies();
    const std::vector<double> & true_beam_traj_Z
        = event.GetTrueTrajZ();
    const std::vector<double> & true_beam_traj_KE
        = event.GetTrueTrajKE();
    const std::vector<int> & true_beam_slices
        = event.GetTrueSlices();
    double beam_inst_P = event.GetBeamInstP();
    const std::vector<double> calibrated_dQdX
        = event.GetdQdXCalibrated();
    const std::vector<double> beam_EField
        = event.GetEField();
    const std::vector<double> track_pitch
        = event.GetTrackPitch();
    const std::vector<double> daughter_thetas
        = event.GetRecoDaughterTrackThetas();

    double beam_shift_delta = 0., beam_energy_delta = 0.;
    if (fBeamShiftCovRoutineActive) {
      beam_shift_delta = GetBeamShiftDelta(true_beam_incidentEnergies);

      double e1 = sqrt(
          std::pow((beam_inst_P + beam_shift_delta)*1.e3, 2) + 139.57*139.57);
      double e2 = sqrt(
          std::pow(beam_inst_P*1.e3, 2) + 139.57*139.57);
      beam_energy_delta = e1 - e2;
      //std::cout << "delta: " << beam_shift_delta << " " << beam_energy_delta <<
      //             std::endl;

      beam_inst_P += beam_shift_delta;
    }


    double end_energy = true_beam_interactingEnergy;
    if (fSliceMethod == "Traj") {
      end_energy = sqrt(true_beam_endP*true_beam_endP*1.e6 + 139.57*139.57) - 139.57;
    }
    else if (fSliceMethod == "E") {
      end_energy = sqrt(true_beam_endP*true_beam_endP*1.e6 + 139.57*139.57) - 139.57;
    }
    else if (fSliceMethod == "Alt") {
      int bin = fEndSlices->GetXaxis()->FindBin(true_beam_traj_Z.back());
      if (bin > 0)
        end_energy = fMeans.at(bin);
    }

    if (samples.find(sample_ID) == samples.end()) {
      std::cout << "Warning: skipping sample " << sample_ID << std::endl;
      continue;
    }

    //Build the true incident energy vector based on the slice cut
    std::vector<double> good_true_incEnergies;
    if (fill_incident) {
      good_true_incEnergies = MakeTrueIncidentEnergies(
        true_beam_traj_Z, true_beam_traj_KE, true_beam_slices,
        true_beam_incidentEnergies);
    }

    int bin = GetBeamBin(beam_energy_bins, (use_beam_inst_P ? beam_inst_P :
                                            true_beam_startP),
                                           (fBeamShiftCovRoutineActive ||
                                            fRestrictBeamInstP));

    if ((fBeamShiftCovRoutineActive ||
         fRestrictBeamInstP)
        && bin == -1) {
      continue;
    }
    //int bin = GetBeamBin(beam_energy_bins, (use_beam_inst_P ? beam_inst_P :
    //                                        true_beam_startP));

    //Weight for the event
    //Possibly affected by signal, flux, and syst parameters
    double weight = 1.;

    //Base class
    std::vector<ThinSliceSample> & samples_vec = samples.at(sample_ID)[bin];
    bool is_signal = signal_sample_checks.at(sample_ID);

    int signal_index = -1;
    ThinSliceSample * this_sample = 0x0;
    if (!is_signal) {
      this_sample = &samples_vec.at(0);
    }
    else {
      //Iterate through the true bins and find the correct one
      bool found = false;

      if (samples_vec.size() == 1) {
        this_sample = &(samples_vec[0]);

        signal_index = 0;

        found = true;

        //If fitting under/overflow: Need to consider here
        //int index = j - 1 + (fit_under_over ? 1 : 0);
        weight *= signal_pars.at(sample_ID)[signal_index];
      }
      else {
        for (size_t j = 1; j < samples_vec.size()-1; ++j) {
          ThinSliceSample & sample = samples_vec.at(j);
          if (sample.CheckInSignalRange(end_energy)) {
            this_sample = &sample;

            //Need to fix this -- will be offset for fit_under_over
            signal_index = j;

            found = true;

            //If fitting under/overflow: Need to consider here
            int index = j - 1 + (fit_under_over ? 1 : 0);
            weight *= signal_pars.at(sample_ID)[index];

            break;
          }
        }
      }
      if (!found) {
        //over/underflow here
        if (end_energy <= samples_vec[1].RangeLowEnd()) {
          this_sample = &samples_vec[0];
          signal_index = 0;
          if (fit_under_over) weight *= signal_pars.at(sample_ID)[0];
          //If fitting under/overflow: Need to consider here
        }
        else if (end_energy >
                 samples_vec[samples_vec.size()-2].RangeHighEnd()) {
          //If fitting under/overflow: Need to consider here
          this_sample = &samples_vec.back();
          if (fit_under_over && !tie_under_over) {
            weight *= signal_pars.at(sample_ID).back();
          }
          else if (fit_under_over && tie_under_over) {
            weight *= signal_pars.at(sample_ID)[0];
          }
        }
        else {
          std::cout << "Warning: could not find true bin " <<
                       end_energy << std::endl;
        }
      }
    }

    int flux_type = this_sample->GetFluxType();
    if (flux_pars.find(flux_type) != flux_pars.end()) {
      weight *= flux_pars.at(flux_type);
    }

    double val[1] = {0};
    
    int new_selection = selection_ID;
    TH1D * selected_hist
        = (TH1D*)this_sample->GetSelectionHists().at(new_selection);
    //if (new_selection == 4) {
    if (std::find(fEndZSelections.begin(), fEndZSelections.end(), new_selection)
        != fEndZSelections.end()) {
      if (selected_hist->FindBin(reco_beam_endZ) == 0) {
        val[0] = selected_hist->GetBinCenter(1);
      }
      else if (selected_hist->FindBin(reco_beam_endZ) >
               selected_hist->GetNbinsX()) {
        val[0] = selected_hist->GetBinCenter(selected_hist->GetNbinsX());
      }
      else {
        val[0] = reco_beam_endZ;
      }
    }
    //else if (new_selection > 4) {
    else if (
        std::find(fOneBinSelections.begin(), fOneBinSelections.end(), new_selection)
        != fOneBinSelections.end()) {
      val[0] = .5;
    }
    else if (reco_beam_incidentEnergies.size()) {

      double energy[1] = {0.};
      if (fFakeResolution < 0) {
        energy[0] = {reco_beam_interactingEnergy + beam_energy_delta};
        if (fDoEnergyFix) {
          for (size_t k = 1; k < reco_beam_incidentEnergies.size(); ++k) {
            double deltaE = ((reco_beam_incidentEnergies)[k-1] -
                             (reco_beam_incidentEnergies)[k]);
            if (fVaryCalibrationFakeData && fFillFakeDataInMain) {
              energy[0] += deltaE;
              if (deltaE*fDataCalibrationFactor < fEnergyFix) {
                energy[0] -= deltaE*fDataCalibrationFactor;
              }
            }
            else {
              if (deltaE > fEnergyFix) {
                energy[0] += deltaE; 
              }
            }
          }
        }
      }
      else {
        //std::cout << "Using fake res: " << fFakeResolution << std::endl;
        if (fUseStoredEnergies) {
          energy[0] = fStoredEnergies[worker_id][energy_index]; 
          //std::cout << "Used stored energy " << energy[0] << std::endl;
          ++energy_index;
        }
        else {
          energy[0] = fRNG.Gaus(end_energy, fFakeResolution);
          fStoredEnergies[worker_id].push_back(energy[0]);
          //std::cout << "Storing " << energy[0] << std::endl;
        }
      }

      if (selected_hist->FindBin(energy[0]) == 0) {
        val[0] = selected_hist->GetBinCenter(1);
      }
      else if (selected_hist->FindBin(energy[0]) > selected_hist->GetNbinsX()) {
        val[0] = selected_hist->GetBinCenter(selected_hist->GetNbinsX());
      }
      else {
        val[0] = energy[0];
      }
    }
    else {
      val[0] = selected_hist->GetBinCenter(1);
    }

    //std::cout << "Weights: " << weight << std::endl;
    //if (fSystematics) {
    if (syst_pars.size()) {
      weight *= fSystematics->GetEventWeight(event, signal_index, syst_pars);
      //std::cout << "\tAfter systs: " << weight << std::endl;
    }

    //if (fG4RWPars) {
    if (g4rw_pars.size()) {
      //weight *= fG4RWPars->GetEventWeight(event, signal_index, g4rw_pars);
      weight *= fSystematics->GetSignalWeight_G4RWCoeffNoPar(event, signal_index);
      weight *= fSystematics->GetSignalWeight_TiedG4RWCoeffNoPar(event, signal_index);
      //std::cout << "\tAfter g4rw: " << weight << std::endl;
    }

    if (fix_factors != 0x0) {
      int bin = fix_factors->at(new_selection)->FindBin(val[0]);
      weight *= fix_factors->at(new_selection)->GetBinContent(bin);
      //std::cout << "\tAfter fix: " << weight << std::endl;
    }

    if (fFillFakeDataInMain) {
      weight *= FakeDataWeight(event);
      //std::cout << "\tAfter fake: " << weight << std::endl;
    }

    //HERE ADD THE VARIED STATS weight
    //std::cout << "Stat Var: " << fStatVar << std::endl;
    if (fStatVar) {
      //std::cout << "statvar " << this_sample->GetStatVarWeight(new_selection, val[0]) << std::endl;
      weight *= this_sample->GetStatVarWeight(new_selection, val[0]);
    }

    std::lock_guard<std::mutex> guard(fRefillMutex);
    this_sample->FillSelectionHist(new_selection, val, weight);
    //Fill the total incident hist with truth info
    if (fill_incident) {
      //this_sample->FillTrueIncidentHist(good_true_incEnergies, weight);
      this_sample->AddIncidentEnergies(good_true_incEnergies, weight);
      if (true_beam_incidentEnergies.size() > 0) {
      this_sample->AddESliceEnergies(
          {(true_beam_incidentEnergies)[0],
           end_energy}, weight);
      }
    }
    this_sample->AddVariedFlux(weight);

    FillExtraHistsMC(*this_sample, event, weight);
  }
}



void protoana::AbsCexDriver::RefillMCSamples(
    const std::vector<ThinSliceEvent> & events,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    std::vector<double> & beam_energy_bins,
    const std::map<int, std::vector<double>> & signal_pars,
    const std::map<int, double> & flux_pars,
    const std::map<std::string, ThinSliceSystematic> & syst_pars,
    const std::map<std::string, ThinSliceSystematic> & g4rw_pars,
    bool fit_under_over, bool tie_under_over, bool use_beam_inst_P,
    bool fill_incident, std::map<int, TH1 *> * fix_factors) {

  CALLGRIND_TOGGLE_COLLECT;

  //Reset all samples
  //Base class
  for (auto it = samples.begin(); it != samples.end(); ++it) {
    for (size_t i = 0; i < it->second.size(); ++i) {
      for (size_t j = 0; j < it->second[i].size(); ++j) {
        it->second[i][j].Reset();
      }
    }
  }

  //Calculate events per workers
  size_t nevents = events.size() / fNWorkers;
  size_t remainder = events.size() % fNWorkers;
  std::vector<size_t> events_split(fNWorkers, nevents);
  while (remainder > fNWorkers) {
    nevents = nevents/fNWorkers;
    remainder = nevents % fNWorkers;

    for (size_t i = 0; i < fNWorkers; ++i) {
      events_split[i] += nevents;
    }
  }
  for (size_t i = 0; i < remainder; ++i) {
    events_split[i] += 1;
  }
  //for (auto evs : events_split) {
  //  std::cout << evs << std::endl;
  //}

  std::vector<std::thread> workers;
  for (auto workerid : ROOT::TSeqI(fNWorkers)) {
    std::thread worker(&AbsCexDriver::RefillSampleLoop, this,
        std::ref(events), std::ref(samples), std::ref(signal_sample_checks), std::ref(beam_energy_bins),
        std::ref(signal_pars), std::ref(flux_pars), std::ref(syst_pars), std::ref(g4rw_pars), fit_under_over,
        tie_under_over, use_beam_inst_P, fill_incident, fix_factors,
        workerid, events_split);
    workers.emplace_back(std::move(worker));
    //workers.emplace_back(&AbsCexDriver::RefillSampleLoop, this,
       // events, samples, signal_sample_checks, beam_energy_bins,
       // signal_pars, flux_pars, syst_pars, fit_under_over,
       // tie_under_over, use_beam_inst_P, fill_incident, fix_factors,
       // workerid, events_split);
  }

  for (auto &&worker : workers) { worker.join();}
  if (fFakeResolution > 0. && !fUseStoredEnergies) fUseStoredEnergies = true; //turn this on so that every time after the first, the same energy is used
  //Threading ideas:
  //see https://root.cern/doc/master/mt304__fillHistos_8C.html
  //and https://root.cern/doc/master/mt201__parallelHistoFill_8C.html
  //
  //  -- split up events into chunks
  //  -- elsewhere define member function to pass as a std::thread to pool of workers
  //  -- That function will handle iterating through events
  //  -- Works mostly as normal until sample->FillSelectionHist, etc.
  //  -- Need to make sure using threadsafe versions of hists within there -- same for AddVariedFlux
  //  -- Maybe separate that out and just have a running sum at the end?

  /*
  for (size_t i = 0; i < events.size(); ++i) {
    const ThinSliceEvent & event = events.at(i);
    int sample_ID = event.GetSampleID();
    int selection_ID = event.GetSelectionID();

    double true_beam_interactingEnergy = event.GetTrueInteractingEnergy();
    double reco_beam_interactingEnergy = event.GetRecoInteractingEnergy();
    double true_beam_endP = event.GetTrueEndP();
    double reco_beam_endZ = event.GetRecoEndZ();
    double true_beam_startP = event.GetTrueStartP();

    const std::vector<double> & reco_beam_incidentEnergies
        = event.GetRecoIncidentEnergies();
    const std::vector<double> & true_beam_incidentEnergies
        = event.GetTrueIncidentEnergies();
    const std::vector<double> & true_beam_traj_Z
        = event.GetTrueTrajZ();
    const std::vector<double> & true_beam_traj_KE
        = event.GetTrueTrajKE();
    const std::vector<int> & true_beam_slices
        = event.GetTrueSlices();
    double beam_inst_P = event.GetBeamInstP();
    const std::vector<double> calibrated_dQdX
        = event.GetdQdXCalibrated();
    const std::vector<double> beam_EField
        = event.GetEField();
    const std::vector<double> track_pitch
        = event.GetTrackPitch();
    const std::vector<double> daughter_thetas
        = event.GetRecoDaughterTrackThetas();

    double beam_shift_delta = 0.;
    if (fBeamShiftCovRoutineActive) {
      beam_shift_delta = GetBeamShiftDelta(true_beam_incidentEnergies);
      beam_inst_P += beam_shift_delta;
    }


    double end_energy = true_beam_interactingEnergy;
    if (fSliceMethod == "Traj") {
      end_energy = sqrt(true_beam_endP*true_beam_endP*1.e6 + 139.57*139.57) - 139.57;
    }
    else if (fSliceMethod == "E") {
      end_energy = sqrt(true_beam_endP*true_beam_endP*1.e6 + 139.57*139.57) - 139.57;
    }
    else if (fSliceMethod == "Alt") {
      int bin = fEndSlices->GetXaxis()->FindBin(true_beam_traj_Z.back());
      if (bin > 0)
        end_energy = fMeans.at(bin);
    }

    if (samples.find(sample_ID) == samples.end()) {
      std::cout << "Warning: skipping sample " << sample_ID << std::endl;
      continue;
    }

    //Build the true incident energy vector based on the slice cut
    std::vector<double> good_true_incEnergies;
    if (fill_incident) {
      good_true_incEnergies = MakeTrueIncidentEnergies(
        true_beam_traj_Z, true_beam_traj_KE, true_beam_slices,
        true_beam_incidentEnergies);
    }

    int bin = GetBeamBin(beam_energy_bins, (use_beam_inst_P ? beam_inst_P :
                                            true_beam_startP),
                                           fBeamShiftCovRoutineActive);

    if (fBeamShiftCovRoutineActive && bin == -1) continue;
    //int bin = GetBeamBin(beam_energy_bins, (use_beam_inst_P ? beam_inst_P :
    //                                        true_beam_startP));

    //Weight for the event
    //Possibly affected by signal, flux, and syst parameters
    double weight = 1.;

    //Base class
    std::vector<ThinSliceSample> & samples_vec = samples.at(sample_ID)[bin];
    bool is_signal = signal_sample_checks.at(sample_ID);

    int signal_index = -1;
    ThinSliceSample * this_sample = 0x0;
    if (!is_signal) {
      this_sample = &samples_vec.at(0);
    }
    else {
      //Iterate through the true bins and find the correct one
      bool found = false;
      for (size_t j = 1; j < samples_vec.size()-1; ++j) {
        ThinSliceSample & sample = samples_vec.at(j);
        if (sample.CheckInSignalRange(end_energy)) {
          this_sample = &sample;

          signal_index = j;

          found = true;

          //If fitting under/overflow: Need to consider here
          int index = j - 1 + (fit_under_over ? 1 : 0);
          weight *= signal_pars.at(sample_ID)[index];

          break;
        }
      }
      if (!found) {
        //over/underflow here
        if (end_energy <= samples_vec[1].RangeLowEnd()) {
          this_sample = &samples_vec[0];
          signal_index = 0;
          if (fit_under_over) weight *= signal_pars.at(sample_ID)[0];
          //If fitting under/overflow: Need to consider here
        }
        else if (end_energy >
                 samples_vec[samples_vec.size()-2].RangeHighEnd()) {
          //If fitting under/overflow: Need to consider here
          this_sample = &samples_vec.back();
          if (fit_under_over && !tie_under_over) {
            weight *= signal_pars.at(sample_ID).back();
          }
          else if (fit_under_over && tie_under_over) {
            weight *= signal_pars.at(sample_ID)[0];
          }
        }
        else {
          std::cout << "Warning: could not find true bin " <<
                       end_energy << std::endl;
        }
      }
    }

    int flux_type = this_sample->GetFluxType();
    if (flux_pars.find(flux_type) != flux_pars.end()) {
      weight *= flux_pars.at(flux_type);
    }

    double val[1] = {0};
    
    int new_selection = selection_ID;
    TH1D * selected_hist
        = (TH1D*)this_sample->GetSelectionHists().at(new_selection);
    if (std::find(fEndZSelections.begin(), fEndZSelections.end(), new_selection)
        != fEndZSelections.end()) {
      if (selected_hist->FindBin(reco_beam_endZ) == 0) {
        val[0] = selected_hist->GetBinCenter(1);
      }
      else if (selected_hist->FindBin(reco_beam_endZ) >
               selected_hist->GetNbinsX()) {
        val[0] = selected_hist->GetBinCenter(selected_hist->GetNbinsX());
      }
      else {
        val[0] = reco_beam_endZ;
      }
    }
    else if (
        std::find(fOneBinSelections.begin(), fOneBinSelections.end(), new_selection)
        != fOneBinSelections.end()) {
      val[0] = .5;
    }
    else if (reco_beam_incidentEnergies.size()) {

      double energy[1] = {0.};
      energy[0] = {reco_beam_interactingEnergy + 1.e3*beam_shift_delta};
      if (fDoEnergyFix) {
        for (size_t k = 1; k < reco_beam_incidentEnergies.size(); ++k) {
          double deltaE = ((reco_beam_incidentEnergies)[k-1] -
                           (reco_beam_incidentEnergies)[k]);
          if (deltaE > fEnergyFix) {
            energy[0] += deltaE; 
          }
        }
      }

      if (selected_hist->FindBin(energy[0]) == 0) {
        val[0] = selected_hist->GetBinCenter(1);
      }
      else if (selected_hist->FindBin(energy[0]) > selected_hist->GetNbinsX()) {
        val[0] = selected_hist->GetBinCenter(selected_hist->GetNbinsX());
      }
      else {
        val[0] = energy[0];
      }
    }
    else {
      val[0] = selected_hist->GetBinCenter(1);
    }

    weight *= fSystematics->GetEventWeight(event, signal_index, syst_pars);

    if (fix_factors != 0x0) {
      int bin = fix_factors->at(new_selection)->FindBin(val[0]);
      weight *= fix_factors->at(new_selection)->GetBinContent(bin);
    }

    this_sample->FillSelectionHist(new_selection, val, weight);

    //Fill the total incident hist with truth info
    if (fill_incident) {
      //this_sample->FillTrueIncidentHist(good_true_incEnergies, weight);
      this_sample->AddIncidentEnergies(good_true_incEnergies, weight);
      if (true_beam_incidentEnergies.size() > 0) {
      this_sample->AddESliceEnergies(
          {(true_beam_incidentEnergies)[0],
           end_energy}, weight);
      }
    }

    this_sample->AddVariedFlux(weight);
  }*/
  CALLGRIND_TOGGLE_COLLECT;
}

void protoana::AbsCexDriver::SetupSysts(
    const std::vector<ThinSliceEvent> & events,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    std::vector<double> & beam_energy_bins,
    const std::map<std::string, ThinSliceSystematic> & pars,
    const std::map<std::string, ThinSliceSystematic> & g4rw_pars,
    TFile & output_file) {
  fSystematics = new PDSPSystematics(events, samples, signal_sample_checks,
                                     beam_energy_bins, pars, g4rw_pars, output_file,
                                     (!fInclusive ? 4 : 2), //Upstream
                                     (!fInclusive ? 6 : 4), (!fInclusive ? 7 : 5), //NoTrack, Decay 
                                     (!fInclusive ? 6 : 4), (!fInclusive? 5 : 3), //Past FV, BeamCut
                                     (!fInclusive ? 4 : 2));//Past FV Selection ID
 
  /*fG4RWPars = new PDSPSystematics(events, samples, signal_sample_checks,
                                  beam_energy_bins, g4rw_pars, output_file,
                                  (!fInclusive ? 4 : 2), //Upstream
                                  (!fInclusive ? 6 : 4), (!fInclusive ? 7 : 5), //NoTrack, Decay 
                                  (!fInclusive ? 6 : 4), (!fInclusive? 5 : 3), //Past FV, BeamCut
                                  (!fInclusive ? 4 : 2));//Past FV Selection ID*/
}

/*
void protoana::AbsCexDriver::SetupSyst_BoxBeam(
    const std::map<std::string, ThinSliceSystematic> & pars) {
  if (pars.find("box_beam_weight") == pars.end()) {
    return;
  }
  fBoxBeamRegions
      = pars.at("box_beam_weight").GetOption<
          std::vector<std::pair<double, double>>>("Regions");
  fBoxBeamFraction = pars.at("box_beam_weight").GetOption<double>("Fraction");
}*/


/*
void protoana::AbsCexDriver::SetupSyst_EffVarWeight(
    const std::map<std::string, ThinSliceSystematic> & pars) {
  if (pars.find("eff_var_weight") == pars.end()) {
    return;
  }
  fEffVarF = pars.at("eff_var_weight").GetOption<double>("F");
  fEffVarCut = pars.at("eff_var_weight").GetOption<double>("Cut");
}*/

/*
void protoana::AbsCexDriver::SetupSyst_LowP(
    const std::map<std::string, ThinSliceSystematic> & pars) {
  if (pars.find("low_p_weight") == pars.end()) {
    return;
  }
  fLowPFractions
      = pars.at("low_p_weight").GetOption<std::vector<double>>("Fractions");
}

void protoana::AbsCexDriver::SetupSyst_NPi0(
    const std::map<std::string, ThinSliceSystematic> & pars) {
  if (pars.find("npi0_weight") == pars.end()) {
    return;
  }
  fNPi0Fractions
      = pars.at("npi0_weight").GetOption<std::vector<double>>("Fractions");
}*/

/*
double protoana::AbsCexDriver::GetSystWeight_EffVar(
    const ThinSliceEvent & event,
    const std::map<std::string, ThinSliceSystematic> & pars) {
  
  if (pars.find("eff_var_weight") == pars.end()) return 1.;

  const int selection_ID = event.GetSelectionID();
  if (selection_ID > 3) return 1.;

  const std::vector<double> daughter_thetas 
      = event.GetRecoDaughterTrackThetas();
  const std::vector<double> track_scores
      = event.GetRecoDaughterTrackScores();
  size_t n = 0;
  for (size_t i = 0; i < track_scores.size(); ++i) {
    if ((track_scores[i] > fEffVarCut) &&
        (daughter_thetas[i]*180./TMath::Pi() < 20.) &&
        (daughter_thetas[i] > -999.))
      ++n;
  }
  //std::cout << selection_ID << " Has " << n << " tracks under 20. degrees. Weight: ";
 
  double weight = 1.;
  if (n > 0) {
    if (selection_ID == 3) {
      weight = (1 - std::pow(fEffVarF*pars.at("eff_var_weight").GetValue(), n))/
               (1 - std::pow(fEffVarF, n));
    }
    else {
      weight = (std::pow(pars.at("eff_var_weight").GetValue(), n));
    }
  }
  //std::cout << weight << std::endl;
  //std::cout << "\tF: " << fEffVarF << " Val: " <<
  //             pars.at("eff_var_weight").GetValue() << std::endl;
  return weight;
}*/

/*
void protoana::AbsCexDriver::SetupSyst_EDivWeight(
    const std::map<std::string, ThinSliceSystematic> & pars) {
  if (pars.find("ediv_weight") == pars.end()) {
    return;
  }
  fEDivF = pars.at("ediv_weight").GetOption<double>("F");
  fEDivCut = pars.at("ediv_weight").GetOption<double>("Cut");
}

void protoana::AbsCexDriver::SetupSyst_NoTrackWeight(
    const std::map<std::string, ThinSliceSystematic> & pars) {
  if (pars.find("no_track_weight") == pars.end()) {
    return;
  }
  fNoTrackF = pars.at("no_track_weight").GetOption<double>("F");
}

void protoana::AbsCexDriver::SetupSyst_BeamEffsWeight(
    const std::map<std::string, ThinSliceSystematic> & pars) {
  std::cout << "Beam effs" << std::endl;
  if (pars.find("beam_cut_weight") != pars.end()) {
    fBeamCutF = pars.at("beam_cut_weight").GetOption<double>("F");
  }
  if (pars.find("no_track_weight") != pars.end()) {
    fNoTrackF = pars.at("no_track_weight").GetOption<double>("F");
  }
}

void protoana::AbsCexDriver::SetupSyst_EndZNoTrackWeight(
    const std::map<std::string, ThinSliceSystematic> & pars) {
  if (pars.find("end_z_no_track_weight") != pars.end()) {
    fEndZNoTrackCut = pars.at("end_z_no_track_weight").GetOption<double>("Cut");

    std::vector<std::pair<int, std::vector<double>>> temp
        = pars.at("end_z_no_track_weight").GetOption
            <std::vector<std::pair<int, std::vector<double>>>>("Fractions");
    fEndZFractions = std::map<int, std::vector<double>>(temp.begin(), temp.end());
    //for (auto it = fEndZFractions.begin(); it != fEndZFractions.end(); ++it) {
    //  std::cout << it->first << " ";
    //  for (auto & f : it->second) std::cout << f << " ";
    //  std::cout << std::endl;
    //}

  }
}*/

/*
double protoana::AbsCexDriver::GetSystWeight_EndZNoTrack(
    const ThinSliceEvent & event,
    int signal_index,
    const std::map<std::string, ThinSliceSystematic> & pars) {
  if (pars.find("end_z_no_track_weight") == pars.end())
    return 1.;

  //Upstream Interactions
  if (event.GetSampleID() == 4) return 1.;

  if (fEndZFractions[event.GetSampleID()].size() == 0) return 1.;

  if (event.GetTrueEndZ() > fEndZNoTrackCut) return 1.;

  double variation = pars.at("end_z_no_track_weight").GetValue();

  double fraction = (signal_index > -1 ?
                     fEndZFractions[event.GetSampleID()][signal_index] :
                     fEndZFractions[event.GetSampleID()].back());

  double val = (event.GetSelectionID() == 6 ? variation :
                (1. - variation*fraction)/(1. - fraction));
  if (val < 0.) std::cout << "endznotrack < 0: " << val << std::endl;
  return val;

  //return pars.at("end_z_no_track_weight").GetValue();
}

double protoana::AbsCexDriver::GetSystWeight_UpstreamInt(
    const ThinSliceEvent & event,
    const std::map<std::string, ThinSliceSystematic> & pars) {
  if (pars.find("upstream_int_weight") == pars.end())
    return 1.;

  return ((event.GetSampleID() == 4) ?
          pars.at("upstream_int_weight").GetValue() :
          1.);
}


void protoana::AbsCexDriver::SetupSyst_BeamMatch(
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

double protoana::AbsCexDriver::GetSystWeight_BeamMatch(
    const ThinSliceEvent & event,
    const std::map<std::string, ThinSliceSystematic> & pars) {
  if (pars.find("beam_match_weight") == pars.end())
    return 1.;

  //Upstream Interactions
  if (event.GetSampleID() == 4) return 1.;

  //No reco track
  if (event.GetSelectionID() == 6) return 1.;

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
  return (matched ? variation : (1. - variation*fraction)/(1. - fraction));
}

double protoana::AbsCexDriver::GetSystWeight_BoxBeam(
    const ThinSliceEvent & event,
    const std::map<std::string, ThinSliceSystematic> & pars) {
  if (pars.find("box_beam_weight") == pars.end()) return 1.;

  //if (event.GetSelectionID() != 5) return 1.;

  double startY = event.GetRecoStartY();

  bool near_box_beam = false;
  for (const auto & region : fBoxBeamRegions) {
    if (region.first < startY && startY < region.second) {
      near_box_beam = true;
      break;
    }
  }

  double variation = pars.at("box_beam_weight").GetValue();
  //return (near_box_beam ? variation :
  //        (1. - variation*fBoxBeamFraction)/(1. - fBoxBeamFraction));
  return (near_box_beam ? variation : 1.);
}

double protoana::AbsCexDriver::GetSystWeight_LowP(
    const ThinSliceEvent & event,
    int signal_index,
    const std::map<std::string, ThinSliceSystematic> & pars) {
  if (pars.find("low_p_weight") == pars.end()) return 1.;

  if (event.GetSampleID() != 3) return 1.;

  //std::vector<double> fractions
  //    = pars.at("low_p_weight").GetOption<std::vector<double>>("Fractions");

  std::vector<int> true_beam_daughter_PDG = event.GetTrueDaughterPDGs();
  std::vector<double> true_beam_daughter_startP
      = event.GetTrueDaughterStartPs();

  bool sub_threshold_pi = true;
  for (size_t j = 0; j < true_beam_daughter_PDG.size(); ++j) {
    if (abs(true_beam_daughter_PDG.at(j)) == 211 &&
        true_beam_daughter_startP.at(j) > .150) {
      sub_threshold_pi = false; 
      break;
    }
  }

  double variation = pars.at("low_p_weight").GetValue();
  double fraction = (signal_index > -1 ?
                     fLowPFractions[signal_index] : fLowPFractions.back());

  return (sub_threshold_pi ? variation :
          (1. - variation*fraction)/(1. - fraction));
}

double protoana::AbsCexDriver::GetSystWeight_NPi0(
    const ThinSliceEvent & event,
    int signal_index,
    const std::map<std::string, ThinSliceSystematic> & pars) {
  if (pars.find("npi0_weight") == pars.end()) return 1.;

  if (event.GetSampleID() != 3) return 1.;

  std::vector<int> true_beam_daughter_PDG = event.GetTrueDaughterPDGs();

  bool npi0 = true;
  int pi0_count = 0;
  for (size_t j = 0; j < true_beam_daughter_PDG.size(); ++j) {
    if (true_beam_daughter_PDG.at(j) == 111) {
      ++pi0_count; 
    }
    else if (abs(true_beam_daughter_PDG.at(j)) == 211) {
      npi0 = false;
      break;
    }
  }

  npi0 = npi0 && (pi0_count > 1);

  double variation = pars.at("npi0_weight").GetValue();
  double fraction = (signal_index > -1 ?
                     fNPi0Fractions[signal_index] : fNPi0Fractions.back());

  return (npi0 ? variation :
          (1. - variation*fraction)/(1. - fraction));
}

double protoana::AbsCexDriver::GetSystWeight_EDiv(
    const ThinSliceEvent & event,
    const std::map<std::string, ThinSliceSystematic> & pars) {
  
  if (pars.find("ediv_weight") == pars.end()) return 1.;

  const int selection_ID = event.GetSelectionID();
  if (selection_ID != 4) return 1.;

  const double endZ = event.GetRecoEndZ();
 
  double weight = 1.;
  double var = pars.at("ediv_weight").GetValue();
  if (endZ < fEDivCut) {
    weight = var;
  }
  else {
    weight = (1. - var*fEDivF)/(1. - fEDivF);
  }

  return weight;
}

double protoana::AbsCexDriver::GetSystWeight_NoTrack(
    const ThinSliceEvent & event,
    const std::map<std::string, ThinSliceSystematic> & pars) {
  
  if (pars.find("no_track_weight") == pars.end()) {
    return 1.;
  }

  const int selection_ID = event.GetSelectionID();

  double weight = 1.;
  double var = pars.at("no_track_weight").GetValue();
  if (selection_ID == 6) {
    weight = var;
  }
  else {
    weight = (1. - var*fNoTrackF)/(1. - fNoTrackF);
  }
  //std::cout << selection_ID << " " << fNoTrackF << " " << var << " " << weight << std::endl;

  return weight;
}

double protoana::AbsCexDriver::GetSystWeight_BeamEffs(
    const ThinSliceEvent & event,
    const std::map<std::string, ThinSliceSystematic> & pars) {
  
  if (pars.find("beam_cut_weight") == pars.end() &&
      pars.find("no_track_weight") == pars.end()) {
    return 1.;
  }
  else if (pars.find("beam_cut_weight") != pars.end() &&
           pars.find("no_track_weight") == pars.end()) {
    const int selection_ID = event.GetSelectionID();

    double weight = 1.;
    double var = pars.at("beam_cut_weight").GetValue();
    if (selection_ID == 5) {
      weight = var;
    }
    else {
      weight = (1. - var*fBeamCutF)/(1. - fBeamCutF);
    }

    return weight;
  }
  else if (pars.find("no_track_weight") != pars.end() &&
           pars.find("beam_cut_weight") == pars.end()) {
    const int selection_ID = event.GetSelectionID();

    double weight = 1.;
    double var = pars.at("no_track_weight").GetValue();
    if (selection_ID == 6) {
      weight = var;
    }
    else {
      weight = (1. - var*fNoTrackF)/(1. - fNoTrackF);
    }

    return weight;
  }
  else {
    const int selection_ID = event.GetSelectionID();

    double weight = 1.;
    double var_no_track = pars.at("no_track_weight").GetValue();
    double var_beam_cut = pars.at("beam_cut_weight").GetValue();
    if (selection_ID == 6) {
      weight = var_no_track;
    }
    else if (selection_ID == 5) {
      weight = var_beam_cut;
    }
    else {
      weight = (1. - var_no_track*fNoTrackF - var_beam_cut*fBeamCutF)/
               (1. - fNoTrackF - fBeamCutF);
    }

    if (weight < 0.) {
      std::cout << "Negative beam cut weight: " << var_no_track << " " <<
                fNoTrackF << " " << var_beam_cut << " " << fBeamCutF <<
                std::endl;
    }
    return weight;
  }
}*/

double protoana::AbsCexDriver::GetFakeWeight_G4RWCoeff(
    const ThinSliceEvent & event/*,
    const std::vector<std::string> & branches,
    const std::vector<double> & vars*/) {
  double weight = 1.;

  for (size_t i = 0; i < fG4RWVars.size(); ++i) {
    weight *= event.GetG4RWCoeffWeight(
        fG4RWBranches[i], fG4RWVars[i]);
  }
  return weight;
}

/*
double protoana::AbsCexDriver::GetSystWeight_G4RWCoeff(
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
  return weight;
}

double protoana::AbsCexDriver::GetSystWeight_G4RW(
    const ThinSliceEvent & event,
    const std::map<std::string, ThinSliceSystematic> & pars,
    const ThinSliceSample & sample,
    int selection_ID, double val) {
  double weight = 1.;
  for (std::string & s : fActiveG4RWSysts) {
    if (std::isnan(pars.at(s).GetValue())) {
      std::string message = "protoana::AbsCexDriver::GetSystWeight_G4RW ";
      message += s;
      message += " has nan value";
      throw std::runtime_error(message);
    }
    weight *= sample.GetSplineWeight(s, pars.at(s).GetValue(), selection_ID, val);
    if (weight < 0.) {
      std::cout << "G4RW turned negative for " << s << std::endl;
    }
  }

  if (weight < 0.) {
    std::cout << "Warning: returning negative value for event " <<
                 event.GetEventID() << " " << event.GetSubrunID() << " " <<
                 event.GetRunID() << std::endl;
  }
  return weight;
}*/

//double protoana::AbsCexDriver::GetSystWeight_BeamShift(
//    const ThinSliceEvent & event,
//    const std::map<std::string, ThinSliceSystematic> & pars) {
//  if (pars.find("beam_shift") == pars.end()) return 1.;
//  if (event.GetPDG() != 211) return 1.;
//  double x_val = pars.at("beam_shift").GetValue();
//  double y_val = (event.GetBeamInstP() - event.GetTrueStartP())/
//                  event.GetTrueStartP();
//  if (y_val < fSystBeamShiftLimits.first/*fSystBeamShiftMap->GetYmin()*/ ||
//      y_val > fSystBeamShiftLimits.second/*fSystBeamShiftMap->GetYmax()*/) {
//    return 1.;
//  }
//
//  double nominal_mean = fSystBeamShiftMeans->Eval(0.);
//  double nominal_width = fSystBeamShiftWidths->Eval(0.);
//  double varied_mean = fSystBeamShiftMeans->Eval(x_val);
//  double varied_width = fSystBeamShiftWidths->Eval(x_val);
//
//  //double weight = 1.;
//
//  /*
//  if (y_val > fSystBeamShiftRatioLimitUp) {
//    if (x_val != 0.) std::cout << x_val << " Past limit: " << y_val << " ";
//    weight = ROOT::Math::normal_cdf_c(y_val, varied_width, varied_mean);
//    if (x_val != 0.) std::cout << weight << " ";
//    weight /= ROOT::Math::normal_cdf_c(y_val, nominal_width, nominal_mean);
//    if (x_val != 0.) std::cout << weight << std::endl;
//  }
//  else if (y_val < fSystBeamShiftRatioLimitDown) {
//    if (x_val != 0.) std::cout << x_val << " Past limit: " << y_val << " ";
//    weight = ROOT::Math::normal_cdf(y_val, varied_width, varied_mean);
//    if (x_val != 0.) std::cout << weight << " ";
//    weight /= ROOT::Math::normal_cdf(y_val, nominal_width, nominal_mean);
//    if (x_val != 0.) std::cout << weight << std::endl;
//  }
//  else {*/
//  double weight = (nominal_width/varied_width)*
//                  exp(.5*std::pow(((y_val - nominal_mean)/nominal_width), 2)
//                      - .5*std::pow(((y_val - varied_mean)/varied_width), 2));
//
//  if (weight > fSystBeamShiftWeightCap && fSystBeamShiftWeightCap > 0.) {
//    //std::cout << "Weight above cap: " << weight << " " << x_val <<
//    //             " event: " << event.GetEventID() << std::endl;
//    weight = fSystBeamShiftWeightCap;
//  }
//  //}
//
//  //fSystBeamShiftWeight = fSystBeamShiftMap->Interpolate(x_val, y_val);
//  if (fSystBeamShiftTreeSave) {
//    fSystBeamShiftWeight = weight;
//    fSystBeamShiftVal = x_val;
//    fSystBeamShiftR = y_val;
//    fSystBeamShiftTree->Fill();
//  }
//  return weight;
//}

/*
double protoana::AbsCexDriver::GetSystWeight_BeamShift2D(
    const ThinSliceEvent & event,
    const std::map<std::string, ThinSliceSystematic> & pars) {
  bool has_shift = pars.find("beam_2D_shift") != pars.end();
  bool has_B = pars.find("beam_2D_B") != pars.end();
  if (!has_shift && !has_B) return 1.;
  if (event.GetPDG() != 211) return 1.;

  double x_val = (has_B ? pars.at("beam_2D_B").GetValue() : 0.);
  double y_val = (has_shift ? pars.at("beam_2D_shift").GetValue() : 0.);

  double r_val = (event.GetBeamInstP() - event.GetTrueStartP())/
                  event.GetTrueStartP();
  if (abs(r_val) > .2) {
    return 1.;
  }
  
  double nominal_mean = fSystBeam2DMeans->Interpolate(0., 0.);
  double nominal_std_dev = fSystBeam2DStdDevs->Interpolate(0., 0.);
  double varied_mean = fSystBeam2DMeans->Interpolate(x_val, y_val);
  double varied_std_dev = fSystBeam2DStdDevs->Interpolate(x_val, y_val);

  double weight = (nominal_std_dev/varied_std_dev)*
                  exp(.5*std::pow(((r_val - nominal_mean)/nominal_std_dev), 2)
                      - .5*std::pow(((r_val - varied_mean)/varied_std_dev), 2));
  if (isnan(weight)) {
    std::cout << "WARNING: RETURNED NAN " << weight << " " <<
                 varied_mean << " " << varied_std_dev << " " << x_val <<
                 " " << y_val << std::endl;
  }

  //if (fSystBeamShiftMap->Interpolate(x_val, y_val) < 0.) {
  //  std::cout << "WARNING: Interpolated " << x_val << " " << y_val <<
  //               fSystBeamShiftMap->Interpolate(x_val, y_val) << std::endl;
  //}
  //else if (isnan(fSystBeamShiftMap->Interpolate(x_val, y_val))) {
  //  std::cout << "WARNING: Interpolated " << x_val << " " << y_val <<
  //               fSystBeamShiftMap->Interpolate(x_val, y_val) << std::endl;
  //}

  fSystBeamShift2DWeight = weight;
  fSystBeamShift2DBVal = x_val;
  fSystBeamShift2DVal = y_val;
  fSystBeamShift2DR = r_val;
  fSystBeamShift2DTree->Fill();
  return weight;//fSystBeamShiftMap->Interpolate(x_val, y_val);
}*/


void protoana::AbsCexDriver::SetupExtraHists(
    ThinSliceDataSet & data_set, 
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & fake_samples) {

  data_set.SetupExtraHists(fExtraHistSets);

  /*for (auto & hist_set : fExtraHistSets) {
    std::string category = hist_set.get<std::string>("Category");
    fActiveExtraHistsData.push_back({category, fStoredExtraHistsData[category]});
    fActiveExtraHistsMC.push_back({category, fStoredExtraHistsMC[category]});
  }*/

  for (auto & hist_set : fExtraHistSets) {
    std::string category = hist_set.get<std::string>("Category");
    if (category == "EndZ") {
      fActiveExtraHistsData.push_back({category, FillExtraHistDataEndZ});
      fActiveExtraHistsMC.push_back({category, FillExtraHistMCEndZ});
    }
    if (category == "EndZGoodReco") {
      fActiveExtraHistsData.push_back({category, FillExtraHistDataEndZGoodReco});
      fActiveExtraHistsMC.push_back({category, FillExtraHistMCEndZGoodReco});
    }
    if (category == "StartX") {
      fActiveExtraHistsData.push_back({category, FillExtraHistDataStartX});
      fActiveExtraHistsMC.push_back({category, FillExtraHistMCStartX});
      auto options = hist_set.get<fhicl::ParameterSet>("Options");
      fStartXDataMean = options.get<double>("DataMean", 0.);
      fStartXMCMean = options.get<double>("MCMean", 0.);
      fStartXDataSigma = options.get<double>("DataSigma", 0.);
      fStartXMCSigma = options.get<double>("MCSigma", 0.);
    }
    if (category == "StartY") {
      fActiveExtraHistsData.push_back({category, FillExtraHistDataStartY});
      fActiveExtraHistsMC.push_back({category, FillExtraHistMCStartY});
      auto options = hist_set.get<fhicl::ParameterSet>("Options");
      fStartYDataMean = options.get<double>("DataMean", 0.);
      fStartYMCMean = options.get<double>("MCMean", 0.);
      fStartYDataSigma = options.get<double>("DataSigma", 0.);
      fStartYMCSigma = options.get<double>("MCSigma", 0.);
    }
    if (category == "StartZ") {
      fActiveExtraHistsData.push_back({category, FillExtraHistDataStartZ});
      fActiveExtraHistsMC.push_back({category, FillExtraHistMCStartZ});
      auto options = hist_set.get<fhicl::ParameterSet>("Options");
      fStartZDataMean = options.get<double>("DataMean", 0.);
      fStartZMCMean = options.get<double>("MCMean", 0.);
      fStartZDataSigma = options.get<double>("DataSigma", 0.);
      fStartZMCSigma = options.get<double>("MCSigma", 0.);
    }
    else if (category == "VertexMichel") {
      fActiveExtraHistsData.push_back({category, FillExtraHistDataVertexMichel});
      fActiveExtraHistsMC.push_back({category, FillExtraHistMCVertexMichel});
    }
    else if (category == "TrackScore") {
      fActiveExtraHistsData.push_back({category, FillExtraHistDataTrackScore});
      fActiveExtraHistsMC.push_back({category, FillExtraHistMCTrackScore});
    }
    else if (category == "Trunc_dEdX") {
      fActiveExtraHistsData.push_back({category, FillExtraHistDataTrunc_dEdX});
      fActiveExtraHistsMC.push_back({category, FillExtraHistMCTrunc_dEdX});
    }
    else if (category == "Chi2PerHit") {
      fActiveExtraHistsData.push_back({category, FillExtraHistDataChi2PerHit});
      fActiveExtraHistsMC.push_back({category, FillExtraHistMCChi2PerHit});
    }
  }

  for (auto & it : samples) {
    int i = 0;
    for (auto & sample_vec : it.second) {
      int j = 0;
      for (auto & sample : sample_vec) {
        sample.SetupExtraHists(fExtraHistSets, i, j);
        ++j;
      }
      ++i;
    }
  }
  for (auto & it : fake_samples) {
    int i = 0;
    for (auto & sample_vec : it.second) {
      int j = 0;
      for (auto & sample : sample_vec) {
        sample.SetupExtraHists(fExtraHistSets, i, j);
        ++j;
      }
      ++i;
    }
  }
}

void protoana::AbsCexDriver::FillExtraHistDataEndZ(ThinSliceDataSet & data_set,
                                               const ExtraHistDataVars & vars) {
  std::string category = "EndZ";
  data_set.FillExtraHist(category, vars.reco_endz);
}
void protoana::AbsCexDriver::FillExtraHistDataEndZGoodReco(ThinSliceDataSet & data_set,
                                               const ExtraHistDataVars & vars) {

  //Skip bad-reco events
  if (vars.selection_ID == 5 || vars.selection_ID == 6) return;

  std::string category = "EndZGoodReco";
  data_set.FillExtraHist(category, vars.reco_endz);
}

void protoana::AbsCexDriver::FillExtraHistDataStartX(ThinSliceDataSet & data_set,
                                               const ExtraHistDataVars & vars) {
  std::string category = "StartX";
  double val = (vars.reco_startx_sce - fStartXDataMean)/fStartXDataSigma;
  data_set.FillExtraHist(category, val);
}
void protoana::AbsCexDriver::FillExtraHistDataStartY(ThinSliceDataSet & data_set,
                                               const ExtraHistDataVars & vars) {
  std::string category = "StartY";
  double val = (vars.reco_starty_sce - fStartYDataMean)/fStartYDataSigma;
  data_set.FillExtraHist(category, val);
}
void protoana::AbsCexDriver::FillExtraHistDataStartZ(ThinSliceDataSet & data_set,
                                               const ExtraHistDataVars & vars) {
  std::string category = "StartZ";
  double val = (vars.reco_startz_sce - fStartZDataMean)/fStartZDataSigma;
  data_set.FillExtraHist(category, val);
}

void protoana::AbsCexDriver::FillExtraHistDataTrackScore(ThinSliceDataSet & data_set,
                                               const ExtraHistDataVars & vars) {
  //Skip non-interaction selected events
  if (vars.selection_ID > 3) return;

  std::string category = "TrackScore";
  for (auto & score : vars.reco_product_track_scores)
    data_set.FillExtraHist(category, score);
}

void protoana::AbsCexDriver::FillExtraHistDataVertexMichel(ThinSliceDataSet & data_set,
                                               const ExtraHistDataVars & vars) {
  if (vars.selection_ID > 3 && vars.selection_ID < 7) return;

  std::string category = "VertexMichel";
  if (vars.vertex_nhits == 0) {
    auto * h = (TH1D*)data_set.GetExtraHist(category);
    h->AddBinContent(1);
  }
  else {
    data_set.FillExtraHist(category,
                           vars.vertex_michel_score/vars.vertex_nhits);
  }
}
void protoana::AbsCexDriver::FillExtraHistDataTrunc_dEdX(
    ThinSliceDataSet & data_set, const ExtraHistDataVars & vars) {
  //Skip non-interaction selected events
  if (vars.selection_ID > 3) return;

  std::string category = "Trunc_dEdX";
  for (auto & val : vars.reco_product_truncated_dEdXs)
    data_set.FillExtraHist(category, val);
}

void protoana::AbsCexDriver::FillExtraHistDataChi2PerHit(
    ThinSliceDataSet & data_set, const ExtraHistDataVars & vars) {
  //Skip non-interaction selected events
  if (vars.selection_ID > 3) return;

  std::string category = "Chi2PerHit";
  for (size_t i = 0; i < vars.reco_product_chi2s_per_hit.size(); ++i) {
    double trunc_dedx = vars.reco_product_truncated_dEdXs[i];
    if (((trunc_dedx < 0.5) || (trunc_dedx > 2.8 && trunc_dedx < 3.4)) &&
        (vars.reco_product_chi2s_per_hit[i] > 0)) {
      data_set.FillExtraHist(category, vars.reco_product_chi2s_per_hit[i]);
    }
  }
}

void protoana::AbsCexDriver::FillExtraHistMCEndZ(ThinSliceSample & sample,
                                                 const ThinSliceEvent & event,
                                                 double weight) {
  std::string category = "EndZ";
  sample.FillExtraHist(category, event.GetRecoEndZ(), weight);
}
void protoana::AbsCexDriver::FillExtraHistMCEndZGoodReco(ThinSliceSample & sample,
                                                 const ThinSliceEvent & event,
                                                 double weight) {
  if (event.GetSelectionID() == 5 || event.GetSelectionID() == 6) return;

  std::string category = "EndZGoodReco";
  sample.FillExtraHist(category, event.GetRecoEndZ(), weight);
}

void protoana::AbsCexDriver::FillExtraHistMCStartX(ThinSliceSample & sample,
                                                 const ThinSliceEvent & event,
                                                 double weight) {
  std::string category = "StartX";
  double val = (event.GetRecoStartX_SCE() - fStartXMCMean)/fStartXMCSigma;
  sample.FillExtraHist(category, val, weight);
}
void protoana::AbsCexDriver::FillExtraHistMCStartY(ThinSliceSample & sample,
                                                 const ThinSliceEvent & event,
                                                 double weight) {
  std::string category = "StartY";
  double val = (event.GetRecoStartY_SCE() - fStartYMCMean)/fStartYMCSigma;
  sample.FillExtraHist(category, val, weight);
}
void protoana::AbsCexDriver::FillExtraHistMCStartZ(ThinSliceSample & sample,
                                                 const ThinSliceEvent & event,
                                                 double weight) {
  std::string category = "StartZ";
  double val = (event.GetRecoStartZ_SCE() - fStartZMCMean)/fStartZMCSigma;
  sample.FillExtraHist(category, val, weight);
}

void protoana::AbsCexDriver::FillExtraHistMCTrackScore(ThinSliceSample & sample,
                                                 const ThinSliceEvent & event,
                                                 double weight) {
  if (event.GetSelectionID() > 3) return;

  std::string category = "TrackScore";
  for (auto & score : event.GetRecoDaughterTrackScores())
    sample.FillExtraHist(category, score, weight);
}

void protoana::AbsCexDriver::FillExtraHistMCVertexMichel(ThinSliceSample & sample,
                                                 const ThinSliceEvent & event,
                                                 double weight) {
  if (event.GetSelectionID() > 3 && event.GetSelectionID() < 7) return;

  std::string category = "VertexMichel";
  if (event.GetVertexNHits() == 0) {
    auto * h = (TH1D*)sample.GetExtraHist(category);
    h->AddBinContent(1, weight);
  }
  else {
    sample.FillExtraHist(
        category,
        event.GetVertexMichelScore()/event.GetVertexNHits(),
        weight);
  }
}

void protoana::AbsCexDriver::FillExtraHistMCTrunc_dEdX(ThinSliceSample & sample,
                                                 const ThinSliceEvent & event,
                                                 double weight) {
  if (event.GetSelectionID() > 3) return;

  std::string category = "Trunc_dEdX";
  for (auto & val : event.GetTrunc_dEdXs())
    sample.FillExtraHist(category, val, weight);
}

void protoana::AbsCexDriver::FillExtraHistMCChi2PerHit(ThinSliceSample & sample,
                                                       const ThinSliceEvent & event,
                                                       double weight) {
  if (event.GetSelectionID() > 3) return;

  std::string category = "Chi2PerHit";
  auto & truncated_dEdXs = event.GetTrunc_dEdXs();
  auto & chi2s_per_hit = event.GetChi2sPerHit();
  for (size_t i = 0; i < chi2s_per_hit.size(); ++i) {
    double trunc_dedx = truncated_dEdXs[i];
    if (((trunc_dedx < 0.5) || (trunc_dedx > 2.8 && trunc_dedx < 3.4))
        && (chi2s_per_hit[i] > 0)) {
      sample.FillExtraHist(category, chi2s_per_hit[i], weight);
    }
  }
}

void protoana::AbsCexDriver::FillExtraHistsData(ThinSliceDataSet & data_set,
                                            const ExtraHistDataVars & vars) {
  for (auto & filler : fActiveExtraHistsData) {
    filler.second(data_set, vars);
  }
}

void protoana::AbsCexDriver::FillExtraHistsMC(ThinSliceSample & sample,
                                              const ThinSliceEvent & event,
                                              double weight) {
  for (auto & filler : fActiveExtraHistsMC) {
    //std::cout << "Filling extra " << filler.first << std::endl;
    filler.second(sample, event, weight);
  }
}

void protoana::AbsCexDriver::BuildDataHists(
    TTree * tree, ThinSliceDataSet & data_set, double & flux,
    const std::vector<double> & beam_energy_bins,
    std::vector<double> & beam_fluxes,
    int split_val) {
  int selection_ID; 
  double reco_beam_interactingEnergy, reco_beam_endZ,
         reco_beam_startX, reco_beam_startY, reco_beam_startZ;
  std::vector<double> * reco_beam_incidentEnergies = 0x0;
  if (!fInclusive) {
    tree->SetBranchAddress("selection_ID", &selection_ID);
  }
  else {
    tree->SetBranchAddress("selection_ID_inclusive", &selection_ID);
  }
  tree->SetBranchAddress("reco_beam_interactingEnergy",
                        &reco_beam_interactingEnergy);
  tree->SetBranchAddress("reco_beam_endZ", &reco_beam_endZ);
  tree->SetBranchAddress("reco_beam_calo_startX", &reco_beam_startX);
  tree->SetBranchAddress("reco_beam_calo_startY", &reco_beam_startY);
  tree->SetBranchAddress("reco_beam_calo_startZ", &reco_beam_startZ);
  tree->SetBranchAddress("reco_beam_incidentEnergies",
                        &reco_beam_incidentEnergies);
  double vertex_michel_score;
  int vertex_nhits;
  tree->SetBranchAddress("reco_beam_vertex_michel_score", &vertex_michel_score);
  tree->SetBranchAddress("reco_beam_vertex_nHits", &vertex_nhits);
  double beam_inst_P;
  tree->SetBranchAddress("beam_inst_P", &beam_inst_P);

  std::vector<double> * reco_daughter_track_scores = 0x0;
  tree->SetBranchAddress("reco_daughter_PFP_trackScore_collection",
                         &reco_daughter_track_scores);

  std::vector<double> * reco_daughter_truncated_dEdXs = 0x0,
                      * reco_daughter_chi2s = 0x0;
  tree->SetBranchAddress("reco_daughter_allTrack_truncLibo_dEdX_pos",
                         &reco_daughter_truncated_dEdXs);
  tree->SetBranchAddress("reco_daughter_allTrack_Chi2_proton",
                         &reco_daughter_chi2s);

  std::vector<int> * reco_daughter_chi2_nhits = 0x0;
  tree->SetBranchAddress("reco_daughter_allTrack_Chi2_ndof",
                         &reco_daughter_chi2_nhits);

  //SetupExtraHists(data_set);
  //data_set.SetupExtraHists(fExtraHistSets);

  TH1D & incident_hist = data_set.GetIncidentHist();
  std::map<int, TH1 *> & selected_hists = data_set.GetSelectionHists();

  if (split_val < 0 || split_val > tree->GetEntries()) {
    flux = tree->GetEntries()/* - split_val*/;
  }
  else {
    flux = split_val;
  }

  double pi_mass_sq = 139.57*139.57;
  for (size_t i = 0; i < beam_energy_bins.size()-1; ++i) {beam_fluxes.push_back(0.);}
  int n_skipped = 0;
  for (int i = 0/*split_val*/; i < flux/*tree->GetEntries()*/; ++i) {
    tree->GetEntry(i);

    double beam_inst_P_scaled = fBeamInstPScale*beam_inst_P;
    int beam_bin = GetBeamBin(beam_energy_bins, beam_inst_P_scaled/*beam_inst_P*/,
                              fRestrictBeamInstP);
    if (beam_bin == -1) {
      ++n_skipped;
      if (fDebugRestrictBeamP)
        std::cout << "skipping " << beam_inst_P_scaled << std::endl;
      continue;
    }

    double deltaE_scale
        = sqrt(std::pow(beam_inst_P_scaled*1.e3, 2) + pi_mass_sq) -
               sqrt(std::pow(beam_inst_P*1.e3, 2) + pi_mass_sq);
    //std::cout << "deltaE_scale: " << deltaE_scale << std::endl; 
    beam_fluxes[beam_bin] += 1.;

    ExtraHistDataVars extra_vars(selection_ID, reco_beam_endZ,
                                 reco_beam_startX, reco_beam_startY,
                                 reco_beam_startZ, vertex_michel_score,
                                 vertex_nhits);
    extra_vars.AddRecoProductTrackScores((*reco_daughter_track_scores));
    extra_vars.AddRecoProductTrunc_dEdXs((*reco_daughter_truncated_dEdXs));

    auto chi2_vec = (*reco_daughter_chi2s);
    for (size_t j = 0; j < chi2_vec.size(); ++j) {
      int nhits = (*reco_daughter_chi2_nhits)[j];
      if (nhits > 0) {
        chi2_vec[j] /= nhits;
      }
      else {
        chi2_vec[j] = -999.;
      }
    } 
    extra_vars.AddRecoProductChi2sPerHit(chi2_vec);

    FillExtraHistsData(data_set, extra_vars);

    double val = 0.;
    if (std::find(fEndZSelections.begin(), fEndZSelections.end(), selection_ID)
        != fEndZSelections.end()) {
      if (selected_hists[selection_ID]->FindBin(reco_beam_endZ) == 0) {
        val = selected_hists[selection_ID]->GetBinCenter(1);
      }
      else if (selected_hists[selection_ID]->FindBin(reco_beam_endZ) >
               selected_hists[selection_ID]->GetNbinsX()) {
        val = selected_hists[selection_ID]->GetBinCenter(
            selected_hists[selection_ID]->GetNbinsX());
      }
      else {
        val = reco_beam_endZ;
      }
    }
    else if (
        std::find(fOneBinSelections.begin(), fOneBinSelections.end(),
                  selection_ID)
        != fOneBinSelections.end()) {
      val = .5;
    }
    else if (reco_beam_incidentEnergies->size()) {
      for (size_t j = 0; j < reco_beam_incidentEnergies->size(); ++j) {
        incident_hist.Fill((*reco_beam_incidentEnergies)[j]);
      }
      if (selected_hists.find(selection_ID) != selected_hists.end()) {
        if (selection_ID != 4 && selection_ID != 5 && selection_ID != 6) {
          double energy = reco_beam_interactingEnergy + deltaE_scale;
          if (fDoEnergyFix) {
            for (size_t k = 1; k < reco_beam_incidentEnergies->size(); ++k) {
              double deltaE = ((*reco_beam_incidentEnergies)[k-1] -
                               (*reco_beam_incidentEnergies)[k]);
              if (fVaryDataCalibration) {
                energy += deltaE;
                if (deltaE*fDataCalibrationFactor < fEnergyFix) {
                  energy -= deltaE*fDataCalibrationFactor;
                }
              }
              else {
                if (deltaE > fEnergyFix) {
                  energy += deltaE; 
                }
              }
            }
          }
          if (selected_hists[selection_ID]->FindBin(energy) == 0) {
            val = selected_hists[selection_ID]->GetBinCenter(1);
          }
          else if (selected_hists[selection_ID]->FindBin(energy) >
                   selected_hists[selection_ID]->GetNbinsX()) {
            val = selected_hists[selection_ID]->GetBinCenter(
                selected_hists[selection_ID]->GetNbinsX());
          }
          else {
            val = energy;
          }
        }
      }
    }
    else {
      val = selected_hists[selection_ID]->GetBinCenter(1);
    }

    selected_hists[selection_ID]->Fill(val);
  }
  std::cout << "Data -- Skipped: " << n_skipped << std::endl;
  flux -= n_skipped;

}

void protoana::AbsCexDriver::BuildFakeData(
    TTree * tree,
    const std::vector<ThinSliceEvent> & events,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    ThinSliceDataSet & data_set, double & flux,
    std::map<int, std::vector<double>> & sample_scales,
    std::vector<double> & beam_energy_bins,
    std::vector<double> & beam_fluxes,
    int split_val, bool scale_to_data_beam_p) {
  //std::string routine = fExtraOptions.get<std::string>("FakeDataRoutine");
  if (fFakeDataRoutine == "SampleScales") {
    FakeDataSampleScales(events, samples, signal_sample_checks, data_set, flux,
                         sample_scales, beam_energy_bins, beam_fluxes,
                         split_val, scale_to_data_beam_p);
  }
  else if (fFakeDataRoutine == "G4RWGrid") {
    FakeDataG4RWGrid(events, samples, signal_sample_checks, data_set, flux,
                 sample_scales, beam_energy_bins, beam_fluxes,
                 split_val, scale_to_data_beam_p);
  }
  else if (fFakeDataRoutine == "EffVar") {
    FakeDataEffVar(events, samples, signal_sample_checks, data_set, flux,
                 sample_scales, split_val);
  }
  else if (fFakeDataRoutine == "LowP") {
    FakeDataLowP(events, samples, signal_sample_checks, data_set, flux,
                 sample_scales, split_val);
  }
  else if (fFakeDataRoutine == "BinnedScales") {
    FakeDataBinnedScales(tree, samples, signal_sample_checks, data_set, flux,
                         sample_scales, split_val);
  }
  else if (fFakeDataRoutine == "dEdX") {
    FakeDatadEdX(tree, samples, signal_sample_checks, data_set, flux,
                 sample_scales, split_val);
  }
  else if (fFakeDataRoutine == "PionAngle") {
    FakeDataPionAngle(tree, samples, signal_sample_checks, data_set, flux,
                      sample_scales, split_val);
  }
  else if (fFakeDataRoutine == "AngleVar") {
    FakeDataAngleVar(events/*tree*/, samples, signal_sample_checks, data_set, flux,
                     beam_energy_bins, beam_fluxes,
                     sample_scales, split_val, scale_to_data_beam_p);
  }
  else if (fFakeDataRoutine == "BeamWeight") {
    FakeDataBeamWeight(events, samples, signal_sample_checks, data_set, flux,
                      sample_scales, split_val);
  }
  else if (fFakeDataRoutine == "BeamScale") {
    std::cout << "Beam scale" << std::endl;
    FakeDataBeamScale(events, samples, signal_sample_checks, data_set, flux,
                      beam_energy_bins, beam_fluxes,
                      sample_scales, split_val, scale_to_data_beam_p);
  }
  /*else if (fFakeDataRoutine == "SetResolution") {
    FakeDataSetResolution(events, samples, signal_sample_checks, data_set, flux,
                     beam_energy_bins, beam_fluxes,
                     sample_scales, split_val, scale_to_data_beam_p);
  }*/
}

void protoana::AbsCexDriver::FakeDataSampleScales(
    const std::vector<ThinSliceEvent> & events,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    ThinSliceDataSet & data_set, double & flux,
    std::map<int, std::vector<double>> & sample_scales,
    std::vector<double> & beam_energy_bins,
    std::vector<double> & beam_fluxes,
    int split_val,
    bool norm_to_data_beam_P) {


  //Build the map for fake data scales
  std::vector<std::pair<int, double>> temp_vec
      = fExtraOptions.get<std::vector<std::pair<int, double>>>(
          "FakeDataScales");
  std::map<int, double> fake_data_scales(temp_vec.begin(), temp_vec.end()); 


  TH1D & incident_hist = data_set.GetIncidentHist();
  std::map<int, TH1 *> & selected_hists = data_set.GetSelectionHists();

  double new_flux = 0.;
  flux = events.size();

  std::map<int, std::vector<double>> nominal_samples;
  for (auto it = sample_scales.begin(); it != sample_scales.end(); ++it) {
    nominal_samples[it->first] = std::vector<double>(it->second.size(), 0.);
  }

  for (size_t i = 0; i < beam_energy_bins.size()-1; ++i) {beam_fluxes.push_back(0.);}

  for (size_t i = 0; i < events.size(); ++i) {

    const ThinSliceEvent & event = events.at(i);
    double true_beam_startP = event.GetTrueStartP();
    double true_beam_endP = event.GetTrueEndP();
    double end_energy = event.GetTrueInteractingEnergy();
    const std::vector<double> & true_beam_traj_Z = event.GetTrueTrajZ();
    const std::vector<double> & true_beam_traj_KE = event.GetTrueTrajKE();
    int sample_ID = event.GetSampleID();
    int selection_ID = event.GetSelectionID();
    const std::vector<double> & true_beam_incidentEnergies
        = event.GetTrueIncidentEnergies();
    const std::vector<double> & reco_beam_incidentEnergies
        = event.GetRecoIncidentEnergies();
    double reco_beam_interactingEnergy = event.GetRecoInteractingEnergy();
    double reco_beam_endZ = event.GetRecoEndZ();
    double beam_P = event.GetBeamInstP();
    const std::vector<int> & true_beam_slices = event.GetTrueSlices();

    //double end_energy = true_beam_interactingEnergy;
    if (fSliceMethod == "Traj") {
      end_energy = sqrt(true_beam_endP*true_beam_endP*1.e6 + 139.57*139.57) - 139.57;
    }
    else if (fSliceMethod == "E") {
      end_energy = sqrt(true_beam_endP*true_beam_endP*1.e6 + 139.57*139.57) - 139.57;
    }
    else if (fSliceMethod == "Alt") {
      int bin = fEndSlices->GetXaxis()->FindBin(true_beam_traj_Z./*->*/back());
      if (bin > 0)
        end_energy = fMeans.at(bin);
    }

    //Add under/overflow check here
    if (samples.find(sample_ID) == samples.end())
      continue;

    double scale = (fake_data_scales.find(sample_ID) != fake_data_scales.end() ?
                    fake_data_scales.at(sample_ID) : 1.);

    int bin = GetBeamBin(beam_energy_bins, (norm_to_data_beam_P ? beam_P :
                                            true_beam_startP));

    //If it's signal
    //Determine if the over/underflow bin 
    bool is_signal = signal_sample_checks.at(sample_ID);
    ThinSliceSample * this_sample = 0x0;
    if (is_signal) {
      std::vector<ThinSliceSample> & samples_vec = samples[sample_ID][bin];
      //Get the samples vec from the first beam energy bin
      bool found = false;

      if (samples_vec.size() == 1) {
        found = true;
        sample_scales[sample_ID][0] += scale;
        nominal_samples[sample_ID][0] += 1.;
        this_sample = &(samples_vec[0]);
      }
      else {
        for (size_t j = 1; j < samples_vec.size()-1; ++j) {
          ThinSliceSample & sample = samples_vec.at(j);
          if (sample.CheckInSignalRange(end_energy)) {     
            found = true;
            sample_scales[sample_ID][j] += scale;
            nominal_samples[sample_ID][j] += 1.;
            this_sample = &sample;
            break;
          }
        }
      }
      if (!found) {//If in the under/overflow, just set to 1. 

        scale = 1.;
        if (end_energy < samples_vec[1].RangeLowEnd()) {
          this_sample = &samples_vec[0];
          sample_scales[sample_ID][0] += scale;
          nominal_samples[sample_ID][0] += 1.;
        }
        else if (end_energy >
                 samples_vec[samples_vec.size()-2].RangeHighEnd()) {

          this_sample = &samples_vec.back();
          sample_scales[sample_ID].back() += scale;
          nominal_samples[sample_ID].back() += 1.;
        }
      }
    }
    else {
      sample_scales[sample_ID][0] += scale;
      nominal_samples[sample_ID][0] += 1.;
      this_sample = &samples[sample_ID][bin][0];
    }

    //this_sample
    new_flux += scale; //1 or scaled

    double val = 0.;
    if (selection_ID == 4) {
      if (selected_hists[selection_ID]->FindBin(reco_beam_endZ) == 0) {
        val = selected_hists[selection_ID]->GetBinCenter(1);
      }
      else if (selected_hists[selection_ID]->FindBin(reco_beam_endZ) >
               selected_hists[selection_ID]->GetNbinsX()) {
        val = selected_hists[selection_ID]->GetBinCenter(
            selected_hists[selection_ID]->GetNbinsX());
      }
      else {
        val = reco_beam_endZ;
      }
    }
    else if (selection_ID > 4) {
      val = .5;
    }
    else if (reco_beam_incidentEnergies.size()) {
      for (size_t j = 0; j < reco_beam_incidentEnergies.size(); ++j) {
        incident_hist.Fill(reco_beam_incidentEnergies[j]);
      }
      if (selected_hists.find(selection_ID) != selected_hists.end()) {
        if (selection_ID != 4 && selection_ID != 5 && selection_ID != 6) {
          double energy = reco_beam_interactingEnergy;
          if (fDoEnergyFix) {
            for (size_t k = 1; k < reco_beam_incidentEnergies.size(); ++k) {
              double deltaE = (reco_beam_incidentEnergies[k-1] -
                               reco_beam_incidentEnergies[k]);
              if (deltaE > fEnergyFix) {
                energy += deltaE; 
              }
            }
          }
          if (selected_hists[selection_ID]->FindBin(energy) == 0) {
            val = selected_hists[selection_ID]->GetBinCenter(1);
          }
          else if (selected_hists[selection_ID]->FindBin(energy) >
                   selected_hists[selection_ID]->GetNbinsX()) {
            val = selected_hists[selection_ID]->GetBinCenter(
                selected_hists[selection_ID]->GetNbinsX());
          }
          else {
            val = energy;
          }
        }
      }
    }
    else {
      val = selected_hists[selection_ID]->GetBinCenter(1);
    }

    selected_hists[selection_ID]->Fill(val, scale);
    std::vector<double> good_true_incEnergies = MakeTrueIncidentEnergies(
        true_beam_traj_Z, true_beam_traj_KE, true_beam_slices,
        true_beam_incidentEnergies);
    this_sample->AddVariedFlux(scale);
    this_sample->AddIncidentEnergies(good_true_incEnergies, scale);
    if (norm_to_data_beam_P) {
      beam_fluxes[GetBeamBin(beam_energy_bins, beam_P)] += scale;
    }
  }

  for (auto it = sample_scales.begin(); it != sample_scales.end(); ++it) {
    for (size_t i = 0; i < it->second.size(); ++i) {
      if (it->second[i] > 0.) {
        it->second[i] /= nominal_samples[it->first][i];
      }
      else {
        it->second[i] = 1.;
      }
      it->second[i] *= (flux/new_flux);
    }
  }

  incident_hist.Scale(flux/new_flux);
  for (auto it = selected_hists.begin(); it != selected_hists.end(); ++it) {
    it->second->Scale(flux/new_flux);
  }
}

void protoana::AbsCexDriver::FakeDataBinnedScales(
    TTree * tree,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    ThinSliceDataSet & data_set, double & flux,
    std::map<int, std::vector<double>> & sample_scales, int split_val) {

  //Build the map for fake data scales
  std::vector<std::pair<int, std::vector<double>>> temp_vec
      = fExtraOptions.get<std::vector<std::pair<int, std::vector<double>>>>(
          "FakeDataBinnedScales");
  std::map<int, std::vector<double>> fake_data_scales(temp_vec.begin(), temp_vec.end()); 

  int sample_ID, selection_ID; 
  double true_beam_interactingEnergy, reco_beam_interactingEnergy;
  double true_beam_endP;
  std::vector<double> * reco_beam_incidentEnergies = 0x0,
                      * true_beam_traj_KE = 0x0,
                      * true_beam_traj_Z = 0x0;
  std::vector<int>    * true_beam_slices = 0x0;
  double reco_beam_endZ;
  tree->SetBranchAddress("reco_beam_endZ", &reco_beam_endZ);
  tree->SetBranchAddress("new_interaction_topology", &sample_ID);
  tree->SetBranchAddress("selection_ID", &selection_ID);
  tree->SetBranchAddress("true_beam_interactingEnergy",
                         &true_beam_interactingEnergy);
  tree->SetBranchAddress("true_beam_endP", &true_beam_endP);
  tree->SetBranchAddress("true_beam_traj_KE", &true_beam_traj_KE);
  tree->SetBranchAddress("true_beam_slices", &true_beam_slices);
  tree->SetBranchAddress("true_beam_traj_Z", &true_beam_traj_Z);
  tree->SetBranchAddress("reco_beam_interactingEnergy",
                         &reco_beam_interactingEnergy);
  tree->SetBranchAddress("reco_beam_incidentEnergies",
                         &reco_beam_incidentEnergies);
  std::vector<double> * true_beam_incidentEnergies = 0x0;
  tree->SetBranchAddress("true_beam_incidentEnergies",
                         &true_beam_incidentEnergies);

  TH1D & incident_hist = data_set.GetIncidentHist();
  std::map<int, TH1 *> & selected_hists = data_set.GetSelectionHists();


  double new_flux = 0.;
  flux = tree->GetEntries() - split_val;

  std::map<int, std::vector<double>> nominal_samples;
  for (auto it = sample_scales.begin(); it != sample_scales.end(); ++it) {
    nominal_samples[it->first] = std::vector<double>(it->second.size(), 0.);
  }

  for (int i = /*0*/split_val; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);

    double end_energy = true_beam_interactingEnergy;
    if (fSliceMethod == "Traj") {
      end_energy = sqrt(true_beam_endP*true_beam_endP*1.e6 + 139.57*139.57) - 139.57;
    }
    else if (fSliceMethod == "E") {
      end_energy = sqrt(true_beam_endP*true_beam_endP*1.e6 + 139.57*139.57) - 139.57;
    }
    else if (fSliceMethod == "Alt") {
      int bin = fEndSlices->GetXaxis()->FindBin(true_beam_traj_Z->back());
      if (bin > 0) {
        end_energy = fMeans.at(bin);
      }
    }

    //Add under/overflow check here
    if (samples.find(sample_ID) == samples.end())
      continue;
    double scale = 1.;

    bool is_scaled = (fake_data_scales.find(sample_ID) !=
                      fake_data_scales.end());

    //If it's signal
    //Determine if the over/underflow bin 
    bool is_signal = signal_sample_checks.at(sample_ID);
    ThinSliceSample * this_sample = 0x0;
    if (is_signal) {
      std::vector<ThinSliceSample> & samples_vec = samples[sample_ID][0];
      //Get the samples vec from the first beam energy bin
      bool found = false;
      for (size_t j = 1; j < samples_vec.size()-1; ++j) {
        ThinSliceSample & sample = samples_vec.at(j);
        if (sample.CheckInSignalRange(end_energy)) {     
          found = true;
          if (is_scaled) {
            scale = fake_data_scales.at(sample_ID)[j];
            sample_scales[sample_ID][j] += scale;
            nominal_samples[sample_ID][j] += 1.;
          }
          this_sample = &sample;
          break;
        }
      }
      if (!found) {
        if (end_energy < samples_vec[1].RangeLowEnd()) {
          this_sample = &samples_vec[0];
          if (is_scaled) {
            scale = fake_data_scales.at(sample_ID)[0];
            sample_scales[sample_ID][0] += scale;
            nominal_samples[sample_ID][0] += 1.;
          }
        }
        else if (end_energy >
                 samples_vec[samples_vec.size()-2].RangeHighEnd()) {
          this_sample = &samples_vec.back();
          if (is_scaled) {
            scale = fake_data_scales.at(sample_ID).back();
            sample_scales[sample_ID].back() += scale;
            nominal_samples[sample_ID].back() += 1.;
          }
        }
      }
    }
    if (!is_signal) {
      if (is_scaled) {
        scale = fake_data_scales.at(sample_ID)[0];
        sample_scales[sample_ID][0] += scale;
        nominal_samples[sample_ID][0] += 1.;
      }
      this_sample = &samples[sample_ID][0][0];
    }


    new_flux += scale; //1 or scaled

    double val = 0.;
    if (selection_ID == 4) {
      if (selected_hists[selection_ID]->FindBin(reco_beam_endZ) == 0) {
        val = selected_hists[selection_ID]->GetBinCenter(1);
      }
      else if (selected_hists[selection_ID]->FindBin(reco_beam_endZ) >
               selected_hists[selection_ID]->GetNbinsX()) {
        val = selected_hists[selection_ID]->GetBinCenter(
            selected_hists[selection_ID]->GetNbinsX());
      }
      else {
        val = reco_beam_endZ;
      }
    }
    else if (selection_ID > 4) {
      val = .5;
    }
    else if (reco_beam_incidentEnergies->size()) {
      for (size_t j = 0; j < reco_beam_incidentEnergies->size(); ++j) {
        incident_hist.Fill((*reco_beam_incidentEnergies)[j]);
      }
      if (selected_hists.find(selection_ID) != selected_hists.end()) {
        if (selection_ID != 4 && selection_ID != 5 && selection_ID != 6) {
          double energy = reco_beam_interactingEnergy;
          if (fDoEnergyFix) {
            for (size_t k = 1; k < reco_beam_incidentEnergies->size(); ++k) {
              double deltaE = ((*reco_beam_incidentEnergies)[k-1] -
                               (*reco_beam_incidentEnergies)[k]);
              if (deltaE > fEnergyFix) {
                energy += deltaE; 
              }
            }
          }
          if (selected_hists[selection_ID]->FindBin(energy) == 0) {
            val = selected_hists[selection_ID]->GetBinCenter(1);
          }
          else if (selected_hists[selection_ID]->FindBin(energy) >
                   selected_hists[selection_ID]->GetNbinsX()) {
            val = selected_hists[selection_ID]->GetBinCenter(
                selected_hists[selection_ID]->GetNbinsX());
          }
          else {
            val = energy;
          }
        }
      }
    }
    else {
      val = selected_hists[selection_ID]->GetBinCenter(1);
    }

    selected_hists[selection_ID]->Fill(val, scale);
    std::vector<double> good_true_incEnergies = MakeTrueIncidentEnergies(
        *true_beam_traj_Z, *true_beam_traj_KE, *true_beam_slices,
        *true_beam_incidentEnergies);
    this_sample->AddVariedFlux(scale);
    this_sample->AddIncidentEnergies(good_true_incEnergies, scale);
  }

  for (auto it = sample_scales.begin(); it != sample_scales.end(); ++it) {
    for (size_t i = 0; i < it->second.size(); ++i) {
      if (it->second[i] > 0.) {
        it->second[i] /= nominal_samples[it->first][i];
      }
      else {
        it->second[i] = 1.;
      }
      it->second[i] *= (flux/new_flux);
    }
  }


  incident_hist.Scale(flux/new_flux);
  for (auto it = selected_hists.begin(); it != selected_hists.end(); ++it) {
    it->second->Scale(flux/new_flux);
  }
}


void protoana::AbsCexDriver::SetupFakeDataG4RW() {
  //Build the map for fake data scales
  fhicl::ParameterSet g4rw_options 
      = fExtraOptions.get<fhicl::ParameterSet>("FakeDataG4RWGrid");

  std::string branch = g4rw_options.get<std::string>("Branch");

  std::vector<size_t> g4rw_pos = g4rw_options.get<std::vector<size_t>>("Position");
  fG4RWVars = g4rw_options.get<std::vector<double>>("Var", {});

  //std::vector<std::string> branches;
  if (g4rw_options.get<bool>("SingleBranch")) {
    fG4RWBranches.push_back(branch);
  }
  else {
    for (size_t & i : g4rw_pos) {
      fG4RWBranches.push_back(branch + "_" + std::to_string(i));
    }
  }

}

void protoana::AbsCexDriver::FakeDataG4RWGrid(
    const std::vector<ThinSliceEvent> & events,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    ThinSliceDataSet & data_set, double & flux,
    std::map<int, std::vector<double>> & sample_scales,
    std::vector<double> & beam_energy_bins,
    std::vector<double> & beam_fluxes,
    int split_val,
    bool norm_to_data_beam_P) {

  //Build the map for fake data scales
  fhicl::ParameterSet g4rw_options 
      = fExtraOptions.get<fhicl::ParameterSet>("FakeDataG4RWGrid");

  std::string branch = g4rw_options.get<std::string>("Branch");

  std::vector<size_t> g4rw_pos = g4rw_options.get<std::vector<size_t>>("Position");
  std::vector<size_t> g4rw_shift = g4rw_options.get<std::vector<size_t>>("Shift");
  std::vector<double> g4rw_var = g4rw_options.get<std::vector<double>>("Var", {});
  bool use_coeffs = g4rw_options.get<bool>("UseCoeffs", false);

  std::vector<std::string> branches;
  if (g4rw_options.get<bool>("SingleBranch")) {
    branches.push_back(branch);
  }
  else {
    for (size_t & i : g4rw_pos) {
      branches.push_back(branch + "_" + std::to_string(i));
    }
  }

  for (size_t i = 0; i < beam_energy_bins.size()-1; ++i) {beam_fluxes.push_back(0.);}

  TH1D & incident_hist = data_set.GetIncidentHist();
  std::map<int, TH1 *> & selected_hists = data_set.GetSelectionHists();

  double new_flux = 0.;
  flux = events.size(); //tree->GetEntries() - split_val;
  

  std::map<int, std::vector<double>> nominal_samples;
  for (auto it = sample_scales.begin(); it != sample_scales.end(); ++it) {
    nominal_samples[it->first] = std::vector<double>(it->second.size(), 0.);
  }

  for (size_t i = 0; i < events.size()/*split_val*/; ++i) {
    const ThinSliceEvent & event = events.at(i);

    int sample_ID = event.GetSampleID();
    int selection_ID = event.GetSelectionID(); 
    double true_beam_interactingEnergy = event.GetTrueInteractingEnergy();
    double reco_beam_interactingEnergy = event.GetRecoInteractingEnergy();
    double true_beam_endP = event.GetTrueEndP();
    double true_beam_startP = event.GetTrueStartP();
    const std::vector<double> & reco_beam_incidentEnergies
        = event.GetRecoIncidentEnergies();
    double reco_beam_endZ = event.GetRecoEndZ();
    const std::vector<double> & true_beam_traj_Z = event.GetTrueTrajZ();
    const std::vector<double> & true_beam_traj_KE = event.GetTrueTrajKE();
    const std::vector<int> & true_beam_slices = event.GetTrueSlices();
    const std::vector<double> & true_beam_incidentEnergies
        = event.GetTrueIncidentEnergies();
    double beam_P = event.GetBeamInstP();


    if (samples.find(sample_ID) == samples.end())
      continue;

    double end_energy = true_beam_interactingEnergy;
    if (fSliceMethod == "Traj") {
      end_energy = sqrt(true_beam_endP*true_beam_endP*1.e6 + 139.57*139.57) - 139.57;
    }
    else if (fSliceMethod == "E") {
      end_energy = sqrt(true_beam_endP*true_beam_endP*1.e6 + 139.57*139.57) - 139.57;
    }
    else if (fSliceMethod == "Alt") {
      int bin = fEndSlices->GetXaxis()->FindBin(true_beam_traj_Z.back());
      if (bin > 0) {
        end_energy = fMeans.at(bin);
      }
    }

    double scale = 1.;

    if (event.GetPDG() == 211) {
    if (!use_coeffs) {
      for (size_t j = 0; j < g4rw_shift.size(); ++j) {
        scale *= event.GetG4RWWeight(branches[j], g4rw_shift[j]);
        //std::cout << branches[j] << " " << g4rw_shift[j] << " " <<
        //             event.GetG4RWWeight(branches[j], g4rw_shift[j]) << " " <<
        //             event.GetEventID() << " " << event.GetSubrunID() << " " <<
        //             event.GetRunID() << std::endl;
      }
    }
    else {
      scale *= GetFakeWeight_G4RWCoeff(event/*, branches, g4rw_var*/);
    }
    }

    //int bin = GetBeamBin(beam_energy_bins, true_beam_startP);
    int bin = GetBeamBin(beam_energy_bins, (norm_to_data_beam_P ? beam_P :
                                            true_beam_startP));

    bool is_signal = signal_sample_checks.at(sample_ID);
    ThinSliceSample * this_sample = 0x0;
    if (is_signal) {
      std::vector<ThinSliceSample> & samples_vec = samples[sample_ID][bin];
      //Get the samples vec from the first beam energy bin
      bool found = false;
      if (samples_vec.size() == 1) {
        found = true;
        sample_scales[sample_ID][0] += scale;
        nominal_samples[sample_ID][0] += 1.;
        this_sample = &(samples_vec[0]);
      }
      else {
        for (size_t j = 1; j < samples_vec.size()-1; ++j) {
          ThinSliceSample & sample = samples_vec.at(j);
          if (sample.CheckInSignalRange(end_energy)) {     
            found = true;
            sample_scales[sample_ID][j] += scale;
            nominal_samples[sample_ID][j] += 1.;
            this_sample = &sample;
            break;
          }
        }
      }

      if (!found) {
        if (end_energy < samples_vec[1].RangeLowEnd()) {
          sample_scales[sample_ID][0] += scale;
          nominal_samples[sample_ID][0] += 1.;
          this_sample = &samples_vec[0];
        }
        else if (end_energy >
                 samples_vec[samples_vec.size()-2].RangeHighEnd()) {
          this_sample = &samples_vec.back();
          sample_scales[sample_ID].back() += scale;
          nominal_samples[sample_ID].back() += 1.;
        }
        else {
         std::cout << "Could not find the sample" << std::endl;
        }
      }
    }
    else {
      this_sample = &samples[sample_ID][bin][0];
      sample_scales[sample_ID][0] += scale;
      nominal_samples[sample_ID][0] += 1.;
    }

    new_flux += scale; //1 or scaled
    double val = 0.;
    if (selection_ID == 4) {
      if (selected_hists[selection_ID]->FindBin(reco_beam_endZ) == 0) {
        val = selected_hists[selection_ID]->GetBinCenter(1);
      }
      else if (selected_hists[selection_ID]->FindBin(reco_beam_endZ) >
               selected_hists[selection_ID]->GetNbinsX()) {
        val = selected_hists[selection_ID]->GetBinCenter(
            selected_hists[selection_ID]->GetNbinsX());
      }
      else {
        val = reco_beam_endZ;
      }
    }
    else if (selection_ID > 4) {
      val = .5;
    }
    else if (reco_beam_incidentEnergies.size()) {
      for (size_t j = 0; j < reco_beam_incidentEnergies.size(); ++j) {
        incident_hist.Fill(reco_beam_incidentEnergies[j]);
      }
      if (selected_hists.find(selection_ID) != selected_hists.end()) {
        if (selection_ID != 4 && selection_ID != 5 && selection_ID != 6) {
          double energy = reco_beam_interactingEnergy;
          if (fDoEnergyFix) {
            for (size_t k = 1; k < reco_beam_incidentEnergies.size(); ++k) {
              double deltaE = (reco_beam_incidentEnergies[k-1] -
                               reco_beam_incidentEnergies[k]);
              if (deltaE > fEnergyFix) {
                energy += deltaE; 
              }
            }
          }
          if (selected_hists[selection_ID]->FindBin(energy) == 0) {
            val = selected_hists[selection_ID]->GetBinCenter(1);
          }
          else if (selected_hists[selection_ID]->FindBin(energy) >
                   selected_hists[selection_ID]->GetNbinsX()) {
            val = selected_hists[selection_ID]->GetBinCenter(
                selected_hists[selection_ID]->GetNbinsX());
          }
          else {
            val = energy;
          }
        }
      }
    }
    else {
      val = selected_hists[selection_ID]->GetBinCenter(1);
    }

    selected_hists[selection_ID]->Fill(val, scale);
    this_sample->AddVariedFlux(scale);
    this_sample->FillSelectionHist(selection_ID, val, scale);
    std::vector<double> good_true_incEnergies = MakeTrueIncidentEnergies(
        true_beam_traj_Z, true_beam_traj_KE, true_beam_slices,
        true_beam_incidentEnergies);
    this_sample->AddIncidentEnergies(good_true_incEnergies, scale);
    if (norm_to_data_beam_P) {
      beam_fluxes[GetBeamBin(beam_energy_bins, beam_P)] += scale;
    }
  }

  for (auto it = sample_scales.begin(); it != sample_scales.end(); ++it) {
    for (size_t i = 0; i < it->second.size(); ++i) {
      if (it->second[i] > 0.) {
        it->second[i] /= nominal_samples[it->first][i];
      }
      else {
        it->second[i] = 1.;
      }
      it->second[i] *= (flux/new_flux);
    }
  }

  incident_hist.Scale(flux/new_flux);
  for (auto it = selected_hists.begin(); it != selected_hists.end(); ++it) {
    it->second->Scale(flux/new_flux);
  }

  std::cout << "Flux/newflux " << flux/new_flux << std::endl;
  ScaleSamples(samples, flux/new_flux);
  //for (auto & samples_vec : samples) {
  //  for (size_t i = 0; i < samples_vec.second.size(); ++i) {
  //    for (size_t j = 0; j < samples_vec.second[j].size(); ++j) {
  //      auto & hists = samples_vec.second[i][j].Get;
  //      for (auto & sel_hist : hists) {
  //        sel_hist.second->Scale(flux/new_flux);
  //      }
  //    }
  //  }
  //}

}

void protoana::AbsCexDriver::ScaleSamples(
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    double scale) {
  for (auto & samples_vec : samples) {
    for (size_t i = 0; i < samples_vec.second.size(); ++i) {
      for (size_t j = 0; j < samples_vec.second[i].size(); ++j) {
        auto & hists = samples_vec.second[i][j].GetSelectionHists();
        for (auto & sel_hist : hists) {
          sel_hist.second->Scale(scale);
        }
      }
    }
  } 
}

void protoana::AbsCexDriver::FakeDataPionAngle(
    TTree * tree,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    ThinSliceDataSet & data_set, double & flux,
    std::map<int, std::vector<double>> & sample_scales, int split_val) {

  //Build the map for fake data scales
  fhicl::ParameterSet options 
      = fExtraOptions.get<fhicl::ParameterSet>("FakeDataPionAngle");

  int sample_ID, selection_ID; 
  double true_beam_interactingEnergy, reco_beam_interactingEnergy;
  double true_beam_endP, true_beam_endPx, true_beam_endPy, true_beam_endPz;
  std::vector<double> * true_beam_daughter_startPx = 0x0,
                      * true_beam_daughter_startPy = 0x0,
                      * true_beam_daughter_startPz = 0x0,
                      * true_beam_daughter_startP = 0x0;
  std::vector<double> * reco_beam_incidentEnergies = 0x0;
  double reco_beam_endZ;
  tree->SetBranchAddress("reco_beam_endZ", &reco_beam_endZ);
  tree->SetBranchAddress("new_interaction_topology", &sample_ID);
  tree->SetBranchAddress("selection_ID", &selection_ID);
  tree->SetBranchAddress("true_beam_interactingEnergy",
                         &true_beam_interactingEnergy);
  tree->SetBranchAddress("true_beam_endP", &true_beam_endP);
  tree->SetBranchAddress("true_beam_endPx", &true_beam_endPx);
  tree->SetBranchAddress("true_beam_endPy", &true_beam_endPy);
  tree->SetBranchAddress("true_beam_endPz", &true_beam_endPz);
  tree->SetBranchAddress("true_beam_daughter_startPx", &true_beam_daughter_startPx);
  tree->SetBranchAddress("true_beam_daughter_startPy", &true_beam_daughter_startPy);
  tree->SetBranchAddress("true_beam_daughter_startPz", &true_beam_daughter_startPz);
  tree->SetBranchAddress("true_beam_daughter_startP",  &true_beam_daughter_startP);
  std::vector<int> * true_beam_daughter_PDG = 0x0;
  tree->SetBranchAddress("true_beam_daughter_PDG", &true_beam_daughter_PDG);
  tree->SetBranchAddress("reco_beam_interactingEnergy",
                         &reco_beam_interactingEnergy);
  tree->SetBranchAddress("reco_beam_incidentEnergies",
                         &reco_beam_incidentEnergies);
  std::vector<double> * true_beam_traj_Z = 0x0,
                      * true_beam_traj_KE = 0x0;
  std::vector<int>    * true_beam_slices = 0x0;
  tree->SetBranchAddress("true_beam_traj_Z", &true_beam_traj_Z);
  tree->SetBranchAddress("true_beam_traj_KE", &true_beam_traj_KE);
  tree->SetBranchAddress("true_beam_slices", &true_beam_slices);
  std::vector<double> * true_beam_incidentEnergies = 0x0;
  tree->SetBranchAddress("true_beam_incidentEnergies",
                         &true_beam_incidentEnergies);

  TH1D & incident_hist = data_set.GetIncidentHist();
  std::map<int, TH1 *> & selected_hists = data_set.GetSelectionHists();

  TFile ratio_file(options.get<std::string>("RatioFile").c_str(), "OPEN");
  std::vector<std::string> ratio_names
      = options.get<std::vector<std::string>>("RatioNames");

  std::vector<TH1D *> ratios;
  for (auto n : ratio_names) {
    ratios.push_back((TH1D*)ratio_file.Get(n.c_str()));
  }
  std::vector<double> limits = options.get<std::vector<double>>("Limits");

  double new_flux = 0.;
  flux = tree->GetEntries() - split_val;
  

  std::map<int, std::vector<double>> nominal_samples;
  for (auto it = sample_scales.begin(); it != sample_scales.end(); ++it) {
    nominal_samples[it->first] = std::vector<double>(it->second.size(), 0.);
  }

  for (int i = split_val; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);

    if (samples.find(sample_ID) == samples.end())
      continue;

    double end_energy = true_beam_interactingEnergy;
    if (fSliceMethod == "Traj") {
      end_energy = sqrt(true_beam_endP*true_beam_endP*1.e6 + 139.57*139.57) - 139.57;
    }
    else if (fSliceMethod == "E") {
      end_energy = sqrt(true_beam_endP*true_beam_endP*1.e6 + 139.57*139.57) - 139.57;
    }
    else if (fSliceMethod == "Alt") {
      int bin = fEndSlices->GetXaxis()->FindBin(true_beam_traj_Z->back());
      if (bin > 0) {
        end_energy = fMeans.at(bin);
      }
    }

    double scale = 1.;
    size_t n_piplus = 0, n_piminus = 0;
    for (size_t j = 0; j < true_beam_daughter_PDG->size(); ++j) {
      if (true_beam_daughter_PDG->at(j) == 211) {
        ++n_piplus;
      }
      else if (true_beam_daughter_PDG->at(j) == -211) {
        ++n_piminus;
      }
    }

    if (sample_ID == 3 && n_piplus == 1 && n_piminus == 0) {
      TH1D * h = 0x0;
      //std::cout << "end p: " << true_beam_endP << std::endl;
      if (true_beam_endP >= limits.back()) {
        h = ratios.back();
        //std::cout << limits.back() << " < " << true_beam_endP << std::endl; 
      }
      else {
        for (size_t j = 1; j < limits.size(); ++j) {
          if (true_beam_endP >= limits[j-1] &&
              true_beam_endP < limits[j]) {
            h = ratios[j-1];
            //std::cout << limits[j-1] << " < " << true_beam_endP <<
            //             " < " << limits[j] << std::endl;
            break;
          }
        }
      }
      //std::cout << "h: " << h << std::endl;

      for (size_t j = 0; j < true_beam_daughter_PDG->size(); ++j) {
        if (true_beam_daughter_PDG->at(j) == 211) {
          double costheta = (true_beam_endPx*true_beam_daughter_startPx->at(j) +
                             true_beam_endPy*true_beam_daughter_startPy->at(j) +
                             true_beam_endPz*true_beam_daughter_startPz->at(j))/
                            (true_beam_endP*true_beam_daughter_startP->at(j));
          
          int bin = h->FindBin(costheta);
          scale *= h->GetBinContent(bin);
        }
      }
    }

    bool is_signal = signal_sample_checks.at(sample_ID);
    ThinSliceSample * this_sample = 0x0;
    if (is_signal) {
      std::vector<ThinSliceSample> & samples_vec = samples[sample_ID][0];
      //Get the samples vec from the first beam energy bin
      bool found = false;
      for (size_t j = 1; j < samples_vec.size()-1; ++j) {
        ThinSliceSample & sample = samples_vec.at(j);
        if (sample.CheckInSignalRange(end_energy)) {     
          found = true;
          sample_scales[sample_ID][j] += scale;
          nominal_samples[sample_ID][j] += 1.;
          this_sample = &sample;
          break;
        }
      }
      if (!found) {
        if (end_energy < samples_vec[1].RangeLowEnd()) {
          sample_scales[sample_ID][0] += scale;
          nominal_samples[sample_ID][0] += 1.;
          this_sample = &samples_vec[0];
        }
        else {
          this_sample = &samples_vec.back();
          sample_scales[sample_ID].back() += scale;
          nominal_samples[sample_ID].back() += 1.;
        }
      }
    }
    else {
      this_sample = &samples[sample_ID][0][0];
      sample_scales[sample_ID][0] += scale;
      nominal_samples[sample_ID][0] += 1.;
    }

    new_flux += scale; //1 or scaled
    double val = 0.;
    if (selection_ID == 4) {
      if (selected_hists[selection_ID]->FindBin(reco_beam_endZ) == 0) {
        val = selected_hists[selection_ID]->GetBinCenter(1);
      }
      else if (selected_hists[selection_ID]->FindBin(reco_beam_endZ) >
               selected_hists[selection_ID]->GetNbinsX()) {
        val = selected_hists[selection_ID]->GetBinCenter(
            selected_hists[selection_ID]->GetNbinsX());
      }
      else {
        val = reco_beam_endZ;
      }
    }
    else if (selection_ID > 4) {
      val = .5;
    }
    else if (reco_beam_incidentEnergies->size()) {
      for (size_t j = 0; j < reco_beam_incidentEnergies->size(); ++j) {
        incident_hist.Fill((*reco_beam_incidentEnergies)[j]);
      }
      if (selected_hists.find(selection_ID) != selected_hists.end()) {
        //if (selection_ID != 4 && selection_ID != 5 && selection_ID != 6) {
        if (selection_ID < 4) {
          double energy = reco_beam_interactingEnergy;
          if (fDoEnergyFix) {
            for (size_t k = 1; k < reco_beam_incidentEnergies->size(); ++k) {
              double deltaE = ((*reco_beam_incidentEnergies)[k-1] -
                               (*reco_beam_incidentEnergies)[k]);
              if (deltaE > fEnergyFix) {
                energy += deltaE; 
              }
            }
          }
          if (selected_hists[selection_ID]->FindBin(energy) == 0) {
            val = selected_hists[selection_ID]->GetBinCenter(1);
          }
          else if (selected_hists[selection_ID]->FindBin(energy) >
                   selected_hists[selection_ID]->GetNbinsX()) {
            val = selected_hists[selection_ID]->GetBinCenter(
                selected_hists[selection_ID]->GetNbinsX());
          }
          else {
            val = energy;
          }
        }
      }
    }
    else {
      val = selected_hists[selection_ID]->GetBinCenter(1);
    }

    selected_hists[selection_ID]->Fill(val, scale);
    this_sample->AddVariedFlux(scale);
    std::vector<double> good_true_incEnergies = MakeTrueIncidentEnergies(
        *true_beam_traj_Z, *true_beam_traj_KE, *true_beam_slices,
        *true_beam_incidentEnergies);
    this_sample->AddIncidentEnergies(good_true_incEnergies, scale);
  }

  for (auto it = sample_scales.begin(); it != sample_scales.end(); ++it) {
    for (size_t i = 0; i < it->second.size(); ++i) {
      if (it->second[i] > 0.) {
        it->second[i] /= nominal_samples[it->first][i];
      }
      else {
        it->second[i] = 1.;
      }
      it->second[i] *= (flux/new_flux);
    }
  }

  incident_hist.Scale(flux/new_flux);
  for (auto it = selected_hists.begin(); it != selected_hists.end(); ++it) {
    it->second->Scale(flux/new_flux);
  }
}

void protoana::AbsCexDriver::FakeDataAngleVar(
    const std::vector<ThinSliceEvent> & events,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    ThinSliceDataSet & data_set, double & flux,
    const std::vector<double> & beam_energy_bins,
    std::vector<double> & beam_fluxes,
    std::map<int, std::vector<double>> & sample_scales, int split_val,
    bool norm_to_data_beam_P) {

  //Build the map for fake data scales
  fhicl::ParameterSet options 
      = fExtraOptions.get<fhicl::ParameterSet>("FakeDataAngleVar");

  TH1D & incident_hist = data_set.GetIncidentHist();
  std::map<int, TH1 *> & selected_hists = data_set.GetSelectionHists();

  TFile * ratio_file = TFile::Open(options.get<std::string>("RatioFile").c_str());
  //std::vector<std::string> ratio_names
  //    = options.get<std::vector<std::string>>("RatioNames");
   auto temp_ratio_names
      = options.get<std::vector<std::pair<int, std::vector<std::string>>>>
          ("RatioNames");
  std::map<int, std::vector<std::string>> ratio_names(
      temp_ratio_names.begin(), temp_ratio_names.end());

  for (size_t i = 0; i < beam_energy_bins.size()-1; ++i) {beam_fluxes.push_back(0.);}

  std::map<int, std::vector<TH1D *>> ratios;
  for (int i = 1; i < 4; ++i) {
    ratios[i] = std::vector<TH1D *>();
    for (auto n : ratio_names[i]) {
      std::string name = n + "_" + std::to_string(i);
      ratios[i].push_back((TH1D*)ratio_file->Get(name.c_str()));
      //std::cout << i << " " << name << " " << ratios[i].back() << std::endl;
    }
  }
  auto temp_limits = options.get<std::vector<std::pair<int, std::vector<double>>>>("Limits");
  std::map<int, std::vector<double>> limits(temp_limits.begin(), temp_limits.end());
  

  double new_flux = 0.;
  flux = events.size();
  

  std::map<int, std::vector<double>> nominal_samples;
  for (auto it = sample_scales.begin(); it != sample_scales.end(); ++it) {
    nominal_samples[it->first] = std::vector<double>(it->second.size(), 0.);
  }

  int check_PDG = options.get<int>("PDG");
  std::vector<double> all_scales;
  for (size_t i = 0; i < events.size()/*split_val*/; ++i) {
    const ThinSliceEvent & event = events.at(i);
    int sample_ID = event.GetSampleID();
    int selection_ID = event.GetSelectionID();
    double beam_P = event.GetBeamInstP();
    double true_beam_interactingEnergy = event.GetTrueInteractingEnergy();
    double reco_beam_interactingEnergy = event.GetRecoInteractingEnergy();
    double true_beam_endP = event.GetTrueEndP();
    const std::vector<double> & reco_beam_incidentEnergies
        = event.GetRecoIncidentEnergies();
    double reco_beam_endZ = event.GetRecoEndZ();

    const std::vector<int> & true_beam_daughter_PDG
        = event.GetTrueDaughterPDGs();
    const std::vector<double> & true_beam_traj_Z = event.GetTrueTrajZ();
    const std::vector<double> & true_beam_traj_KE = event.GetTrueTrajKE();
    const std::vector<int> & true_beam_slices = event.GetTrueSlices();
    const std::vector<double> & true_beam_incidentEnergies
        = event.GetTrueIncidentEnergies();
    double leading_costheta = 0.;
    if (check_PDG == 2212) {
      leading_costheta = event.GetLeadingPCostheta();
    }
    else if (check_PDG == 211) {
      leading_costheta = event.GetLeadingPiPlusCostheta();
    }
    else if (check_PDG == 111) {
      leading_costheta = event.GetLeadingPi0Costheta();
    }

    if (samples.find(sample_ID) == samples.end())
      continue;

    double end_energy = true_beam_interactingEnergy;
    if (fSliceMethod == "Traj") {
      end_energy = sqrt(true_beam_endP*true_beam_endP*1.e6 + 139.57*139.57) - 139.57;
    }
    else if (fSliceMethod == "E") {
      end_energy = sqrt(true_beam_endP*true_beam_endP*1.e6 + 139.57*139.57) - 139.57;
    }
    else if (fSliceMethod == "Alt") {
      int bin = fEndSlices->GetXaxis()->FindBin(true_beam_traj_Z.back());
      if (bin > 0) {
        end_energy = fMeans.at(bin);
      }
    }

    double scale = 1.;
    size_t n_piplus = 0, n_pi0 = 0, n_p = 0;
    for (size_t j = 0; j < true_beam_daughter_PDG.size(); ++j) {
      if (true_beam_daughter_PDG.at(j) == 211) {
        ++n_piplus;
      }
      else if (true_beam_daughter_PDG.at(j) == 111) {
        ++n_pi0;
      }
      else if (true_beam_daughter_PDG.at(j) == 2212) {
        ++n_p;
      }
    }

    double end_p = true_beam_endP*1.e3;
    if (sample_ID < 4 &&
        ((check_PDG == 211 && n_piplus > 0) ||
         (check_PDG == 2212 && n_p > 0) ||
         (check_PDG == 111 && n_pi0 > 0))) {
      TH1D * h = 0x0;
      //std::cout << "end p: " << end_p << " sample_ID: " << sample_ID << std::endl;
      if (end_p >= limits[sample_ID].back()) {
        h = ratios[sample_ID].back();
        //std::cout << limits[sample_ID].back() << " < " << end_p << std::endl; 
      }
      else {
        for (size_t j = 1; j < limits[sample_ID].size(); ++j) {
          if (end_p >= limits[sample_ID][j-1] &&
              end_p < limits[sample_ID][j]) {
            h = ratios[sample_ID][j-1];
            //std::cout << limits[sample_ID][j-1] << " < " << end_p <<
            //             " < " << limits[sample_ID][j] << std::endl;
            break;
          }
        }
      }
      //std::cout << "h: " << h << std::endl;

      int bin = h->FindBin(leading_costheta);
      scale *= h->GetBinContent(bin);
      //std::cout << scale << std::endl;
    }
    all_scales.push_back(scale);

    bool is_signal = signal_sample_checks.at(sample_ID);
    int beam_bin = (norm_to_data_beam_P ?
                    GetBeamBin(beam_energy_bins, beam_P) : 0);
    ThinSliceSample * this_sample = 0x0;
    if (is_signal) {
      std::vector<ThinSliceSample> & samples_vec = samples[sample_ID][beam_bin];
      //Get the samples vec from the first beam energy bin
      bool found = false;
      for (size_t j = 1; j < samples_vec.size()-1; ++j) {
        ThinSliceSample & sample = samples_vec.at(j);
        if (sample.CheckInSignalRange(end_energy)) {     
          found = true;
          sample_scales[sample_ID][j] += scale;
          nominal_samples[sample_ID][j] += 1.;
          this_sample = &sample;
          break;
        }
      }
      if (!found) {
        if (end_energy < samples_vec[1].RangeLowEnd()) {
          sample_scales[sample_ID][0] += scale;
          nominal_samples[sample_ID][0] += 1.;
          this_sample = &samples_vec[0];
        }
        else {
          this_sample = &samples_vec.back();
          sample_scales[sample_ID].back() += scale;
          nominal_samples[sample_ID].back() += 1.;
        }
      }
    }
    else {
      this_sample = &samples[sample_ID][beam_bin][0];
      sample_scales[sample_ID][0] += scale;
      nominal_samples[sample_ID][0] += 1.;
    }

    new_flux += scale; //1 or scaled
    double val = 0.;
    if (selection_ID == 4) {
      if (selected_hists[selection_ID]->FindBin(reco_beam_endZ) == 0) {
        val = selected_hists[selection_ID]->GetBinCenter(1);
      }
      else if (selected_hists[selection_ID]->FindBin(reco_beam_endZ) >
               selected_hists[selection_ID]->GetNbinsX()) {
        val = selected_hists[selection_ID]->GetBinCenter(
            selected_hists[selection_ID]->GetNbinsX());
      }
      else {
        val = reco_beam_endZ;
      }
    }
    else if (selection_ID > 4) {
      val = .5;
    }
    else if (reco_beam_incidentEnergies.size()) {
      for (size_t j = 0; j < reco_beam_incidentEnergies.size(); ++j) {
        incident_hist.Fill((reco_beam_incidentEnergies)[j]);
      }
      if (selected_hists.find(selection_ID) != selected_hists.end()) {
        if (selection_ID < 4) {
          double energy = reco_beam_interactingEnergy;
          if (fDoEnergyFix) {
            for (size_t k = 1; k < reco_beam_incidentEnergies.size(); ++k) {
              double deltaE = ((reco_beam_incidentEnergies)[k-1] -
                               (reco_beam_incidentEnergies)[k]);
              if (deltaE > fEnergyFix) {
                energy += deltaE; 
              }
            }
          }
          if (selected_hists[selection_ID]->FindBin(energy) == 0) {
            val = selected_hists[selection_ID]->GetBinCenter(1);
          }
          else if (selected_hists[selection_ID]->FindBin(energy) >
                   selected_hists[selection_ID]->GetNbinsX()) {
            val = selected_hists[selection_ID]->GetBinCenter(
                selected_hists[selection_ID]->GetNbinsX());
          }
          else {
            val = energy;
          }
        }
      }
    }
    else {
      val = selected_hists[selection_ID]->GetBinCenter(1);
    }

    selected_hists[selection_ID]->Fill(val, scale);
    this_sample->FillSelectionHist(selection_ID, val, scale);
    this_sample->AddVariedFlux(scale);
    std::vector<double> good_true_incEnergies = MakeTrueIncidentEnergies(
        true_beam_traj_Z, true_beam_traj_KE, true_beam_slices,
        true_beam_incidentEnergies);
    this_sample->AddIncidentEnergies(good_true_incEnergies, scale);
    if (norm_to_data_beam_P) {
      beam_fluxes[beam_bin] += scale;
    }
  }

  /* mean is set but not used
  double mean = 0.;
  for (auto s : all_scales) mean += s;
  mean /= all_scales.size();
  //std::cout << "Mean scale: " << mean << std::endl;
  */

  for (auto it = sample_scales.begin(); it != sample_scales.end(); ++it) {
    for (size_t i = 0; i < it->second.size(); ++i) {
      if (it->second[i] > 0.) {
        it->second[i] /= nominal_samples[it->first][i];
      }
      else {
        it->second[i] = 1.;
      }
      it->second[i] *= (flux/new_flux);
    }
  }

  incident_hist.Scale(flux/new_flux);
  for (auto it = selected_hists.begin(); it != selected_hists.end(); ++it) {
    it->second->Scale(flux/new_flux);
  }
  ScaleSamples(samples, flux/new_flux);
}


/*void protoana::AbsCexDriver::FakeDataSetResolution(
    const std::vector<ThinSliceEvent> & events,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    ThinSliceDataSet & data_set, double & flux,
    const std::vector<double> & beam_energy_bins,
    std::vector<double> & beam_fluxes,
    std::map<int, std::vector<double>> & sample_scales, int split_val,
    bool norm_to_data_beam_P) {
 
  fhicl::ParameterSet options 
      = fExtraOptions.get<fhicl::ParameterSet>("FakeDataSetResolution");
  double resolution = options.get<double>("Resolution");

  //TH1D & incident_hist = data_set.GetIncidentHist();
  std::map<int, TH1 *> & selected_hists = data_set.GetSelectionHists();

  for (size_t i = 0; i < beam_energy_bins.size()-1; ++i) {beam_fluxes.push_back(0.);}

  std::map<int, std::vector<TH1D *>> ratios;
  for (int i = 1; i < 4; ++i) {
    ratios[i] = std::vector<TH1D *>();
    for (auto n : ratio_names[i]) {
      std::string name = n + "_" + std::to_string(i);
      ratios[i].push_back((TH1D*)ratio_file->Get(name.c_str()));
      //std::cout << i << " " << name << " " << ratios[i].back() << std::endl;
    }
  }

  double new_flux = 0.;
  flux = events.size();
  

  std::map<int, std::vector<double>> nominal_samples;
  for (auto it = sample_scales.begin(); it != sample_scales.end(); ++it) {
    nominal_samples[it->first] = std::vector<double>(it->second.size(), 0.);
  }

  for (size_t i = 0; i < events.size(); ++i) {
    const ThinSliceEvent & event = events.at(i);
    int sample_ID = event.GetSampleID();
    int selection_ID = event.GetSelectionID();
    double beam_P = event.GetBeamInstP();
    double true_beam_interactingEnergy = event.GetTrueInteractingEnergy();
    double reco_beam_interactingEnergy = event.GetRecoInteractingEnergy();
    double true_beam_endP = event.GetTrueEndP();
    const std::vector<double> & reco_beam_incidentEnergies
        = event.GetRecoIncidentEnergies();
    double reco_beam_endZ = event.GetRecoEndZ();

    const std::vector<int> & true_beam_daughter_PDG
        = event.GetTrueDaughterPDGs();
    const std::vector<double> & true_beam_traj_Z = event.GetTrueTrajZ();
    const std::vector<double> & true_beam_traj_KE = event.GetTrueTrajKE();
    const std::vector<int> & true_beam_slices = event.GetTrueSlices();
    const std::vector<double> & true_beam_incidentEnergies
        = event.GetTrueIncidentEnergies();
    double leading_costheta = 0.;
    if (check_PDG == 2212) {
      leading_costheta = event.GetLeadingPCostheta();
    }
    else if (check_PDG == 211) {
      leading_costheta = event.GetLeadingPiPlusCostheta();
    }
    else if (check_PDG == 111) {
      leading_costheta = event.GetLeadingPi0Costheta();
    }

    if (samples.find(sample_ID) == samples.end())
      continue;

    double end_energy = true_beam_interactingEnergy;
    if (fSliceMethod == "Traj") {
      end_energy = sqrt(true_beam_endP*true_beam_endP*1.e6 + 139.57*139.57) - 139.57;
    }
    else if (fSliceMethod == "E") {
      end_energy = sqrt(true_beam_endP*true_beam_endP*1.e6 + 139.57*139.57) - 139.57;
    }
    else if (fSliceMethod == "Alt") {
      int bin = fEndSlices->GetXaxis()->FindBin(true_beam_traj_Z.back());
      if (bin > 0) {
        end_energy = fMeans.at(bin);
      }
    }

    double scale = 1.;

    bool is_signal = signal_sample_checks.at(sample_ID);
    int beam_bin = (norm_to_data_beam_P ?
                    GetBeamBin(beam_energy_bins, beam_P) : 0);
    ThinSliceSample * this_sample = 0x0;
    if (is_signal) {
      std::vector<ThinSliceSample> & samples_vec = samples[sample_ID][beam_bin];
      //Get the samples vec from the first beam energy bin
      bool found = false;
      for (size_t j = 1; j < samples_vec.size()-1; ++j) {
        ThinSliceSample & sample = samples_vec.at(j);
        if (sample.CheckInSignalRange(end_energy)) {     
          found = true;
          sample_scales[sample_ID][j] += scale;
          nominal_samples[sample_ID][j] += 1.;
          this_sample = &sample;
          break;
        }
      }
      if (!found) {
        if (end_energy < samples_vec[1].RangeLowEnd()) {
          sample_scales[sample_ID][0] += scale;
          nominal_samples[sample_ID][0] += 1.;
          this_sample = &samples_vec[0];
        }
        else {
          this_sample = &samples_vec.back();
          sample_scales[sample_ID].back() += scale;
          nominal_samples[sample_ID].back() += 1.;
        }
      }
    }
    else {
      this_sample = &samples[sample_ID][beam_bin][0];
      sample_scales[sample_ID][0] += scale;
      nominal_samples[sample_ID][0] += 1.;
    }

    new_flux += scale; //1 or scaled

    double val = 0.;
    if (selection_ID == 4) {
      if (selected_hists[selection_ID]->FindBin(reco_beam_endZ) == 0) {
        val = selected_hists[selection_ID]->GetBinCenter(1);
      }
      else if (selected_hists[selection_ID]->FindBin(reco_beam_endZ) >
               selected_hists[selection_ID]->GetNbinsX()) {
        val = selected_hists[selection_ID]->GetBinCenter(
            selected_hists[selection_ID]->GetNbinsX());
      }
      else {
        val = reco_beam_endZ;
      }
    }
    else if (selection_ID > 4) {
      val = .5;
    }
    else if (reco_beam_incidentEnergies.size()) {
      //for (size_t j = 0; j < reco_beam_incidentEnergies.size(); ++j) {
      //  incident_hist.Fill((reco_beam_incidentEnergies)[j]);
      //}
      if (selected_hists.find(selection_ID) != selected_hists.end()) {
        if (selection_ID < 4) {
          val = fRNG.Gaus(end_energy, resolution);
        }
      }
    }
    else {
      val = selected_hists[selection_ID]->GetBinCenter(1);
    }

    selected_hists[selection_ID]->Fill(val, scale);
    this_sample->FillSelectionHist(selection_ID, val, scale);
    this_sample->AddVariedFlux(scale);
    std::vector<double> good_true_incEnergies = MakeTrueIncidentEnergies(
        true_beam_traj_Z, true_beam_traj_KE, true_beam_slices,
        true_beam_incidentEnergies);
    this_sample->AddIncidentEnergies(good_true_incEnergies, scale);
    if (norm_to_data_beam_P) {
      beam_fluxes[beam_bin] += scale;
    }
  }

  for (auto it = sample_scales.begin(); it != sample_scales.end(); ++it) {
    for (size_t i = 0; i < it->second.size(); ++i) {
      if (it->second[i] > 0.) {
        it->second[i] /= nominal_samples[it->first][i];
      }
      else {
        it->second[i] = 1.;
      }
      it->second[i] *= (flux/new_flux);
    }
  }

  //incident_hist.Scale(flux/new_flux);
  for (auto it = selected_hists.begin(); it != selected_hists.end(); ++it) {
    it->second->Scale(flux/new_flux);
  }
  ScaleSamples(samples, flux/new_flux);
}*/

void protoana::AbsCexDriver::FakeDataBeamWeight(
    const std::vector<ThinSliceEvent> & events,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    ThinSliceDataSet & data_set, double & flux,
    std::map<int, std::vector<double>> & sample_scales, int split_val) {

  //Build the map for fake data scales
  fhicl::ParameterSet options 
      = fExtraOptions.get<fhicl::ParameterSet>("FakeDataBeamWeight");


  TH1D & incident_hist = data_set.GetIncidentHist();
  std::map<int, TH1 *> & selected_hists = data_set.GetSelectionHists();

  TFile ratio_file(options.get<std::string>("RatioFile").c_str(), "OPEN");
  TH1D * ratio = (TH1D*)ratio_file.Get("r");

  double new_flux = 0.;
  flux = events.size(); //tree->GetEntries() - split_val;

  std::map<int, std::vector<double>> nominal_samples;
  for (auto it = sample_scales.begin(); it != sample_scales.end(); ++it) {
    nominal_samples[it->first] = std::vector<double>(it->second.size(), 0.);
  }

  for (size_t i = 0; i < events.size(); ++i) {
    const ThinSliceEvent & event = events.at(i);
    int sample_ID = event.GetSampleID();
    int selection_ID = event.GetSelectionID(); 
    double true_beam_interactingEnergy = event.GetTrueInteractingEnergy();
    double reco_beam_interactingEnergy = event.GetRecoInteractingEnergy();
    double true_beam_endP = event.GetTrueEndP();
    const std::vector<double> & reco_beam_incidentEnergies
        = event.GetRecoIncidentEnergies();
    double reco_beam_endZ = event.GetRecoEndZ();
    const std::vector<double> & true_beam_traj_Z = event.GetTrueTrajZ();
    const std::vector<double> & true_beam_traj_KE = event.GetTrueTrajKE();
    const std::vector<int> & true_beam_slices = event.GetTrueSlices();
    const std::vector<double> & true_beam_incidentEnergies
        = event.GetTrueIncidentEnergies();
    double beam_inst_P = event.GetBeamInstP();



    if (samples.find(sample_ID) == samples.end())
      continue;

    double end_energy = true_beam_interactingEnergy;
    if (fSliceMethod == "Traj") {
      end_energy = sqrt(true_beam_endP*true_beam_endP*1.e6 + 139.57*139.57) - 139.57;
    }
    else if (fSliceMethod == "E") {
      end_energy = sqrt(true_beam_endP*true_beam_endP*1.e6 + 139.57*139.57) - 139.57;
    }
    else if (fSliceMethod == "Alt") {
      int bin = fEndSlices->GetXaxis()->FindBin(true_beam_traj_Z.back());
      if (bin > 0) {
        end_energy = fMeans.at(bin);
      }
    }

    int bin = ratio->FindBin(beam_inst_P*1.e3);
    double scale = ratio->GetBinContent(bin);;
    bool is_signal = signal_sample_checks.at(sample_ID);
    ThinSliceSample * this_sample = 0x0;
    if (is_signal) {
      std::vector<ThinSliceSample> & samples_vec = samples[sample_ID][0];
      //Get the samples vec from the first beam energy bin
      bool found = false;
      for (size_t j = 1; j < samples_vec.size()-1; ++j) {
        ThinSliceSample & sample = samples_vec.at(j);
        if (sample.CheckInSignalRange(end_energy)) {     
          found = true;
          sample_scales[sample_ID][j] += scale;
          nominal_samples[sample_ID][j] += 1.;
          this_sample = &sample;
          break;
        }
      }
      if (!found) {
        if (end_energy < samples_vec[1].RangeLowEnd()) {
          sample_scales[sample_ID][0] += scale;
          nominal_samples[sample_ID][0] += 1.;
          this_sample = &samples_vec[0];
        }
        else {
          this_sample = &samples_vec.back();
          sample_scales[sample_ID].back() += scale;
          nominal_samples[sample_ID].back() += 1.;
        }
      }
    }
    else {
      this_sample = &samples[sample_ID][0][0];
      sample_scales[sample_ID][0] += scale;
      nominal_samples[sample_ID][0] += 1.;
    }

    new_flux += scale; //1 or scaled
    double val = 0.;
    if (selection_ID == 4) {
      if (selected_hists[selection_ID]->FindBin(reco_beam_endZ) == 0) {
        val = selected_hists[selection_ID]->GetBinCenter(1);
      }
      else if (selected_hists[selection_ID]->FindBin(reco_beam_endZ) >
               selected_hists[selection_ID]->GetNbinsX()) {
        val = selected_hists[selection_ID]->GetBinCenter(
            selected_hists[selection_ID]->GetNbinsX());
      }
      else {
        val = reco_beam_endZ;
      }
    }
    else if (selection_ID > 4) {
      val = .5;
    }
    else if (reco_beam_incidentEnergies.size()) {
      for (size_t j = 0; j < reco_beam_incidentEnergies.size(); ++j) {
        incident_hist.Fill((reco_beam_incidentEnergies)[j]);
      }
      if (selected_hists.find(selection_ID) != selected_hists.end()) {
        //if (selection_ID != 4 && selection_ID != 5 && selection_ID != 6) {
        if (selection_ID < 4) {
          double energy = reco_beam_interactingEnergy;
          if (fDoEnergyFix) {
            for (size_t k = 1; k < reco_beam_incidentEnergies.size(); ++k) {
              double deltaE = ((reco_beam_incidentEnergies)[k-1] -
                               (reco_beam_incidentEnergies)[k]);
              if (deltaE > fEnergyFix) {
                energy += deltaE; 
              }
            }
          }
          if (selected_hists[selection_ID]->FindBin(energy) == 0) {
            val = selected_hists[selection_ID]->GetBinCenter(1);
          }
          else if (selected_hists[selection_ID]->FindBin(energy) >
                   selected_hists[selection_ID]->GetNbinsX()) {
            val = selected_hists[selection_ID]->GetBinCenter(
                selected_hists[selection_ID]->GetNbinsX());
          }
          else {
            val = energy;
          }
        }
      }
    }
    else {
      val = selected_hists[selection_ID]->GetBinCenter(1);
    }

    selected_hists[selection_ID]->Fill(val, scale);
    this_sample->AddVariedFlux(scale);
    std::vector<double> good_true_incEnergies = MakeTrueIncidentEnergies(
        true_beam_traj_Z, true_beam_traj_KE, true_beam_slices,
        true_beam_incidentEnergies);
    this_sample->AddIncidentEnergies(good_true_incEnergies, scale);
  }

  for (auto it = sample_scales.begin(); it != sample_scales.end(); ++it) {
    for (size_t i = 0; i < it->second.size(); ++i) {
      if (it->second[i] > 0.) {
        it->second[i] /= nominal_samples[it->first][i];
      }
      else {
        it->second[i] = 1.;
      }
      it->second[i] *= (flux/new_flux);
    }
  }

  incident_hist.Scale(flux/new_flux);
  for (auto it = selected_hists.begin(); it != selected_hists.end(); ++it) {
    it->second->Scale(flux/new_flux);
  }
}

void protoana::AbsCexDriver::FakeDatadEdX(
    TTree * tree,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    ThinSliceDataSet & data_set, double & flux,
    std::map<int, std::vector<double>> & sample_scales, int split_val) {

  fhicl::ParameterSet dEdX_options 
      = fExtraOptions.get<fhicl::ParameterSet>("FakeDatadEdX");
  fhicl::ParameterSet cal_set
      = dEdX_options.get<fhicl::ParameterSet>("Cal_set");
  double betaP = cal_set.get<double>("betap");
  double rho = cal_set.get<double>("Rho");
  double wion = cal_set.get<double>("Wion");
  double alpha = cal_set.get<double>("alpha");

  std::vector<fhicl::ParameterSet> PlanePars
      = cal_set.get<std::vector<fhicl::ParameterSet>>("PlaneParameters");
  bool found_collection = false;
  double nominal_CCal = 1.;
  for (auto & p : PlanePars) {
    if (p.get<int>("PlaneID") == 2) {
      nominal_CCal = p.get<double>("calib_factor");
      found_collection = true;
      break;
    }
  }
  
  if (!found_collection) {
    std::string message = "Could not find collection plane calibration factor";
    throw std::runtime_error(message);
  }

  //double nominal_CCal = dEdX_options.get<double>("NominalCCal");
  double varied_CCal = dEdX_options.get<double>("VariedCCal");

  int selection_ID, sample_ID;
  double beam_inst_P;
  double reco_beam_endZ;
  std::vector<double> * calibrated_dQdX = 0x0, * beam_EField = 0x0,
                      * track_pitch = 0x0, * reco_beam_incidentEnergies = 0x0;
  tree->SetBranchAddress("selection_ID", &selection_ID);
  tree->SetBranchAddress("new_interaction_topology", &sample_ID);
  tree->SetBranchAddress("reco_beam_calibrated_dQdX_SCE", &calibrated_dQdX);
  tree->SetBranchAddress("reco_beam_EField_SCE", &beam_EField);
  tree->SetBranchAddress("reco_beam_TrkPitch_SCE", &track_pitch);
  tree->SetBranchAddress("reco_beam_incidentEnergies",
                         &reco_beam_incidentEnergies);
  tree->SetBranchAddress("beam_inst_P", &beam_inst_P);
  tree->SetBranchAddress("reco_beam_endZ", &reco_beam_endZ);

  std::vector<double> * true_beam_traj_Z = 0x0,
                      * true_beam_traj_KE = 0x0;
  std::vector<int>    * true_beam_slices = 0x0;
  tree->SetBranchAddress("true_beam_traj_Z", &true_beam_traj_Z);
  tree->SetBranchAddress("true_beam_traj_KE", &true_beam_traj_KE);
  tree->SetBranchAddress("true_beam_slices", &true_beam_slices);
  std::vector<double> * true_beam_incidentEnergies = 0x0;
  tree->SetBranchAddress("true_beam_incidentEnergies",
                         &true_beam_incidentEnergies);
  double true_beam_endP;
  double true_beam_interactingEnergy;
  tree->SetBranchAddress("true_beam_interactingEnergy",
                         &true_beam_interactingEnergy);
  tree->SetBranchAddress("true_beam_endP", &true_beam_endP);

  std::map<int, TH1 *> & selection_hists = data_set.GetSelectionHists();
  flux = tree->GetEntries() - split_val;
  for (int i = /*0*/split_val; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);
    TH1D * selection_hist = (TH1D*)selection_hists[selection_ID];

    double val = 1.;
    if (selection_ID == 4) {
      if (selection_hist->FindBin(reco_beam_endZ) == 0) {
        val = selection_hist->GetBinCenter(1);
      }
      else if (selection_hist->FindBin(reco_beam_endZ) >
               selection_hist->GetNbinsX()) {
        val = selection_hist->GetBinCenter(
            selection_hist->GetNbinsX());
      }
      else {
        val = reco_beam_endZ;
      }
    }
    else if (selection_ID > 4) {
      val = .5;
    }
    else if (reco_beam_incidentEnergies->size()) {
      double energy = sqrt(beam_inst_P*beam_inst_P*1.e6 + 139.57*139.57) -
                      139.57;
      for (size_t k = 0; k < calibrated_dQdX->size()-1; ++k) {
        if ((*calibrated_dQdX)[k] < 0.) continue;

        double dedx = (nominal_CCal/varied_CCal);
        dedx *= (*calibrated_dQdX)[k];
        dedx *= (betaP / ( rho * (*beam_EField)[k] ) * wion);
        dedx = exp(dedx);
        dedx -= alpha;
        dedx *= ((rho*(*beam_EField)[k])/betaP);

        if (dedx*(*track_pitch)[k] > fEnergyFix)
          continue;
        energy -= dedx*(*track_pitch)[k];
        //std::cout << "Energy: " << energy[0] << " dedx: " << dedx <<
        //             std::endl;
      }
      if (selection_hist->FindBin(energy) == 0) {
        val = selection_hist->GetBinCenter(1);
      }
      else if (selection_hist->FindBin(energy) >
               selection_hist->GetNbinsX()) {
        val = selection_hist->GetBinCenter(selection_hist->GetNbinsX());
      }
      else {
        val = energy;
      }

    }
    else {
      val = selection_hist->GetBinCenter(1);
    }

    double end_energy = true_beam_interactingEnergy;
    if (fSliceMethod == "Traj") {
      end_energy = sqrt(true_beam_endP*true_beam_endP*1.e6 + 139.57*139.57) - 139.57;
    }
    selection_hist->Fill(val);
    bool is_signal = signal_sample_checks.at(sample_ID);
    ThinSliceSample * this_sample = 0x0;
    if (is_signal) {
      std::vector<ThinSliceSample> & samples_vec = samples[sample_ID][0];
      //Get the samples vec from the first beam energy bin
      bool found = false;
      for (size_t j = 1; j < samples_vec.size()-1; ++j) {
        ThinSliceSample & sample = samples_vec.at(j);
        if (sample.CheckInSignalRange(end_energy)) {     
          this_sample = &sample;
          break;
        }
      }
      if (!found) {
        if (end_energy < samples_vec[1].RangeLowEnd()) {
          this_sample = &samples_vec[0];
        }
        else if (end_energy >
                 samples_vec[samples_vec.size()-2].RangeHighEnd()) {
          this_sample = &samples_vec.back();
        }
      }
    }
    else {
      this_sample = &samples[sample_ID][0][0];
    }
    this_sample->AddVariedFlux();
    std::vector<double> good_true_incEnergies = MakeTrueIncidentEnergies(
        *true_beam_traj_Z, *true_beam_traj_KE, *true_beam_slices,
        *true_beam_incidentEnergies);
    this_sample->AddIncidentEnergies(good_true_incEnergies);

  }
  for (auto it = sample_scales.begin(); it != sample_scales.end(); ++it) {
    for (double & s : it->second) {
      s = 1.;
    }
  }
}

void protoana::AbsCexDriver::FakeDataEffVar(
    const std::vector<ThinSliceEvent> & events,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    ThinSliceDataSet & data_set, double & flux,
    std::map<int, std::vector<double>> & sample_scales, int split_val) {

  //Build the map for fake data scales
  fhicl::ParameterSet options 
      = fExtraOptions.get<fhicl::ParameterSet>("FakeDataEffVar");


  TH1D & incident_hist = data_set.GetIncidentHist();
  std::map<int, TH1 *> & selected_hists = data_set.GetSelectionHists();

 // double new_flux = 0.;
  flux = events.size();//tree->GetEntries() - split_val;
  

  std::map<int, std::vector<double>> nominal_samples;
  for (auto it = sample_scales.begin(); it != sample_scales.end(); ++it) {
    nominal_samples[it->first] = std::vector<double>(it->second.size(), 0.);
  }

  double check_val = options.get<double>("F");
  for (size_t i = 0; i < events.size(); ++i) {
    const ThinSliceEvent & event = events.at(i);

    int sample_ID = event.GetSampleID();
    int selection_ID = event.GetSelectionID();
    double true_beam_interactingEnergy
        = event.GetTrueInteractingEnergy();
    double reco_beam_interactingEnergy
        = event.GetRecoInteractingEnergy();
    double true_beam_endP = event.GetTrueEndP();
    const std::vector<double> & reco_beam_incidentEnergies
        = event.GetRecoIncidentEnergies();
    double reco_beam_endZ = event.GetRecoEndZ();
    const std::vector<double> & true_beam_traj_Z = event.GetTrueTrajZ();
    const std::vector<double> & true_beam_traj_KE = event.GetTrueTrajKE();
    const std::vector<int>    & true_beam_slices = event.GetTrueSlices();
    const std::vector<double> & true_beam_incidentEnergies
        = event.GetTrueIncidentEnergies();

    const std::vector<double> & daughter_Theta
        = event.GetRecoDaughterTrackThetas();
    //const std::vector<int> & daughter_true_PDG
    //    = event.GetTrueDaughterPDGs();
    const std::vector<double> & daughter_track_score
        = event.GetRecoDaughterTrackThetas();

    if (samples.find(sample_ID) == samples.end())
      continue;

    double end_energy = true_beam_interactingEnergy;
    if (fSliceMethod == "Traj") {
      end_energy = sqrt(true_beam_endP*true_beam_endP*1.e6 + 139.57*139.57) - 139.57;
    }

    double val = 0.;
    if (selection_ID == 4) {
      if (selected_hists[selection_ID]->FindBin(reco_beam_endZ) == 0) {
        val = selected_hists[selection_ID]->GetBinCenter(1);
      }
      else if (selected_hists[selection_ID]->FindBin(reco_beam_endZ) >
               selected_hists[selection_ID]->GetNbinsX()) {
        val = selected_hists[selection_ID]->GetBinCenter(
            selected_hists[selection_ID]->GetNbinsX());
      }
      else {
        val = reco_beam_endZ;
      }
    }
    else if (selection_ID > 4) {
      val = .5;
    }
    else if (reco_beam_incidentEnergies.size()) {
      for (size_t j = 0; j < reco_beam_incidentEnergies.size(); ++j) {
        incident_hist.Fill(reco_beam_incidentEnergies[j]);
      }
      if (selected_hists.find(selection_ID) != selected_hists.end()) {
        if (selection_ID != 4 && selection_ID != 5 && selection_ID != 6) {
          double energy = reco_beam_interactingEnergy;
          if (fDoEnergyFix) {
            for (size_t k = 1; k < reco_beam_incidentEnergies.size(); ++k) {
              double deltaE = (reco_beam_incidentEnergies[k-1] -
                               reco_beam_incidentEnergies[k]);
              if (deltaE > fEnergyFix) {
                energy += deltaE; 
              }
            }
          }
          if (selected_hists[selection_ID]->FindBin(energy) == 0) {
            val = selected_hists[selection_ID]->GetBinCenter(1);
          }
          else if (selected_hists[selection_ID]->FindBin(energy) >
                   selected_hists[selection_ID]->GetNbinsX()) {
            val = selected_hists[selection_ID]->GetBinCenter(
                selected_hists[selection_ID]->GetNbinsX());
          }
          else {
            val = energy;
          }
        }
      }
    }
    else {
      val = selected_hists[selection_ID]->GetBinCenter(1);
    }

    bool is_signal = signal_sample_checks.at(sample_ID);
    ThinSliceSample * this_sample = 0x0;
    if (is_signal) {
      std::vector<ThinSliceSample> & samples_vec = samples[sample_ID][0];
      //Get the samples vec from the first beam energy bin
      bool found = false;
      for (size_t j = 1; j < samples_vec.size()-1; ++j) {
        ThinSliceSample & sample = samples_vec.at(j);
        if (sample.CheckInSignalRange(end_energy)) {     
          this_sample = &sample;
          break;
        }
      }
      if (!found) {
        if (end_energy < samples_vec[1].RangeLowEnd()) {
          this_sample = &samples_vec[0];
        }
        else if (end_energy >
                 samples_vec[samples_vec.size()-2].RangeHighEnd()) {
          this_sample = &samples_vec.back();
        }
      }
    }
    else {
      this_sample = &samples[sample_ID][0][0];
    }
    this_sample->AddVariedFlux();
    std::vector<double> good_true_incEnergies = MakeTrueIncidentEnergies(
        true_beam_traj_Z, true_beam_traj_KE, true_beam_slices,
        true_beam_incidentEnergies);
    this_sample->AddIncidentEnergies(good_true_incEnergies);

    if (selection_ID < 3) {
      int new_selection = selection_ID;
      for (size_t j = 0; j < daughter_Theta.size(); ++j) {
        if (/*(abs((*daughter_true_PDG)[j]) == 211) &&*/
            (daughter_track_score[j] < 1. && daughter_track_score[j] > .3) &&
            (daughter_Theta[j] > -999) &&
            (daughter_Theta[j]*180./TMath::Pi() < 20.)) {
          double r = fRNG.Uniform();
          if (r < check_val) {
            new_selection = 3;
            break;
          }
        }
      }
      selected_hists[new_selection]->Fill(val);
    }
    else {
      selected_hists[selection_ID]->Fill(val);
    }
  }

  for (auto it = sample_scales.begin(); it != sample_scales.end(); ++it) {
    for (size_t i = 0; i < it->second.size(); ++i) {
      it->second[i] = 1.;
    }
  }
}

void protoana::AbsCexDriver::FakeDataLowP(
    const std::vector<ThinSliceEvent> & events,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    ThinSliceDataSet & data_set, double & flux,
    std::map<int, std::vector<double>> & sample_scales, int split_val) {

  //Build the map for fake data scales
  fhicl::ParameterSet options 
      = fExtraOptions.get<fhicl::ParameterSet>("FakeDataLowP");

  std::vector<double> fractions
      = options.get<std::vector<double>>("Fractions");
  double variation = options.get<double>("Scale");

  TH1D & incident_hist = data_set.GetIncidentHist();
  std::map<int, TH1 *> & selected_hists = data_set.GetSelectionHists();

  double new_flux = 0.;
  flux = events.size();
  

  std::map<int, std::vector<double>> nominal_samples;
  for (auto it = sample_scales.begin(); it != sample_scales.end(); ++it) {
    nominal_samples[it->first] = std::vector<double>(it->second.size(), 0.);
  }

  std::vector<double> all_scales;
  for (size_t i = 0; i < events.size(); ++i) {
    const ThinSliceEvent & event = events.at(i);

    int sample_ID = event.GetSampleID();
    int selection_ID = event.GetSelectionID();
    double true_beam_interactingEnergy = event.GetTrueInteractingEnergy();
    double reco_beam_interactingEnergy = event.GetRecoInteractingEnergy();
    double true_beam_endP = event.GetTrueEndP();
    const std::vector<double> & true_beam_daughter_startP
        = event.GetTrueDaughterStartPs();
    const std::vector<double> & reco_beam_incidentEnergies
        = event.GetRecoIncidentEnergies();
    double reco_beam_endZ = event.GetRecoEndZ();
    const std::vector<int> & true_beam_daughter_PDG
        = event.GetTrueDaughterPDGs();
    const std::vector<double> & true_beam_traj_Z = event.GetTrueTrajZ();
    const std::vector<double> & true_beam_traj_KE = event.GetTrueTrajKE();
    const std::vector<int>    & true_beam_slices = event.GetTrueSlices();
    const std::vector<double> & true_beam_incidentEnergies
        = event.GetTrueIncidentEnergies();

    if (samples.find(sample_ID) == samples.end())
      continue;

    double end_energy = true_beam_interactingEnergy;
    if (fSliceMethod == "Traj") {
      end_energy = sqrt(true_beam_endP*true_beam_endP*1.e6 + 139.57*139.57) - 139.57;
    }
    else if (fSliceMethod == "E") {
      end_energy = sqrt(true_beam_endP*true_beam_endP*1.e6 + 139.57*139.57) - 139.57;
    }
    else if (fSliceMethod == "Alt") {
      int bin = fEndSlices->GetXaxis()->FindBin(true_beam_traj_Z.back());
      if (bin > 0) {
        end_energy = fMeans.at(bin);
      }
    }

    double scale = 1.;

    bool sub_threshold_pi = true;
    for (size_t j = 0; j < true_beam_daughter_PDG.size(); ++j) {
      if (abs(true_beam_daughter_PDG.at(j)) == 211 &&
          true_beam_daughter_startP.at(j) > .150) {
        sub_threshold_pi = false; 
        break;
      }
    }

    bool is_signal = signal_sample_checks.at(sample_ID);
    ThinSliceSample * this_sample = 0x0;
    if (is_signal) {
      std::vector<ThinSliceSample> & samples_vec = samples[sample_ID][0];
      //Get the samples vec from the first beam energy bin
      bool found = false;
      for (size_t j = 1; j < samples_vec.size()-1; ++j) {
        ThinSliceSample & sample = samples_vec.at(j);
        if (sample.CheckInSignalRange(end_energy)) {     

          if (sample_ID == 3 && sub_threshold_pi) {
            scale = variation;
          }
          else if (sample_ID == 3 && !sub_threshold_pi) {
            scale = (1. - variation*fractions[j])/(1. - fractions[j]);
          }

          all_scales.push_back(scale);
          found = true;
          sample_scales[sample_ID][j] += scale;
          nominal_samples[sample_ID][j] += 1.;
          this_sample = &sample;
          break;
        }
      }
      if (!found) {
        if (end_energy < samples_vec[1].RangeLowEnd()) {
          if (sample_ID == 3 && sub_threshold_pi) {
            scale = variation;
          }
          else if (sample_ID == 3 && !sub_threshold_pi) {
            scale = (1. - variation*fractions[0])/(1. - fractions[0]);
          }

          all_scales.push_back(scale);
          sample_scales[sample_ID][0] += scale;
          nominal_samples[sample_ID][0] += 1.;
          this_sample = &samples_vec[0];
        }
        else {
          if (sample_ID == 3 && sub_threshold_pi) {
            scale = variation;
          }
          else if (sample_ID == 3 && !sub_threshold_pi) {
            scale = (1. - variation*fractions.back())/(1. - fractions.back());
          }

          all_scales.push_back(scale);
          this_sample = &samples_vec.back();
          sample_scales[sample_ID].back() += scale;
          nominal_samples[sample_ID].back() += 1.;
        }
      }
    }
    else {
      this_sample = &samples[sample_ID][0][0];
      sample_scales[sample_ID][0] += scale;
      nominal_samples[sample_ID][0] += 1.;
    }

    new_flux += scale; //1 or scaled
    double val = 0.;
    if (selection_ID == 4) {
      if (selected_hists[selection_ID]->FindBin(reco_beam_endZ) == 0) {
        val = selected_hists[selection_ID]->GetBinCenter(1);
      }
      else if (selected_hists[selection_ID]->FindBin(reco_beam_endZ) >
               selected_hists[selection_ID]->GetNbinsX()) {
        val = selected_hists[selection_ID]->GetBinCenter(
            selected_hists[selection_ID]->GetNbinsX());
      }
      else {
        val = reco_beam_endZ;
      }
    }
    else if (selection_ID > 4) {
      val = .5;
    }
    else if (reco_beam_incidentEnergies.size()) {
      for (size_t j = 0; j < reco_beam_incidentEnergies.size(); ++j) {
        incident_hist.Fill(reco_beam_incidentEnergies[j]);
      }
      if (selected_hists.find(selection_ID) != selected_hists.end()) {
        //if (selection_ID != 4 && selection_ID != 5 && selection_ID != 6) {
        if (selection_ID < 4) {
          double energy = reco_beam_interactingEnergy;
          if (fDoEnergyFix) {
            for (size_t k = 1; k < reco_beam_incidentEnergies.size(); ++k) {
              double deltaE = (reco_beam_incidentEnergies[k-1] -
                               reco_beam_incidentEnergies[k]);
              if (deltaE > fEnergyFix) {
                energy += deltaE; 
              }
            }
          }
          if (selected_hists[selection_ID]->FindBin(energy) == 0) {
            val = selected_hists[selection_ID]->GetBinCenter(1);
          }
          else if (selected_hists[selection_ID]->FindBin(energy) >
                   selected_hists[selection_ID]->GetNbinsX()) {
            val = selected_hists[selection_ID]->GetBinCenter(
                selected_hists[selection_ID]->GetNbinsX());
          }
          else {
            val = energy;
          }
        }
      }
    }
    else {
      val = selected_hists[selection_ID]->GetBinCenter(1);
    }

    selected_hists[selection_ID]->Fill(val, scale);
    this_sample->AddVariedFlux(scale);
    std::vector<double> good_true_incEnergies = MakeTrueIncidentEnergies(
        true_beam_traj_Z, true_beam_traj_KE, true_beam_slices,
        true_beam_incidentEnergies);
    this_sample->AddIncidentEnergies(good_true_incEnergies, scale);
  }

  double mean = 0.;
  for (auto s : all_scales) mean += s;
  mean /= all_scales.size();
  std::cout << "Mean scale: " << mean << std::endl;

  for (auto it = sample_scales.begin(); it != sample_scales.end(); ++it) {
    for (size_t i = 0; i < it->second.size(); ++i) {
      if (it->second[i] > 0.) {
        it->second[i] /= nominal_samples[it->first][i];
      }
      else {
        it->second[i] = 1.;
      }
      it->second[i] *= (flux/new_flux);
    }
  }

  incident_hist.Scale(flux/new_flux);
  for (auto it = selected_hists.begin(); it != selected_hists.end(); ++it) {
    it->second->Scale(flux/new_flux);
  }
}

void protoana::AbsCexDriver::FakeDataBeamScale(
    const std::vector<ThinSliceEvent> & events,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    ThinSliceDataSet & data_set, double & flux,
    const std::vector<double> & beam_energy_bins,
    std::vector<double> & beam_fluxes,
    std::map<int, std::vector<double>> & sample_scales, int split_val,
    bool norm_to_data_beam_P) {

  //Build the map for fake data scales
  fhicl::ParameterSet options 
      = fExtraOptions.get<fhicl::ParameterSet>("FakeDataBeamScale");

  std::vector<double> variations = options.get<std::vector<double>>("Scales");
  for (const auto & v: variations) std::cout << v << " ";
  std::cout << std::endl;
  std::vector<double> bin_edges = options.get<std::vector<double>>("Bins");
  for (size_t i = 0; i < beam_energy_bins.size()-1; ++i) {beam_fluxes.push_back(0.);}

  TH1D & incident_hist = data_set.GetIncidentHist();
  std::map<int, TH1 *> & selected_hists = data_set.GetSelectionHists();

  double new_flux = 0.;
  flux = events.size();
  

  std::map<int, std::vector<double>> nominal_samples;
  for (auto it = sample_scales.begin(); it != sample_scales.end(); ++it) {
    nominal_samples[it->first] = std::vector<double>(it->second.size(), 0.);
  }

  std::vector<double> all_scales;
  for (size_t i = 0; i < events.size(); ++i) {
    const ThinSliceEvent & event = events.at(i);

    int sample_ID = event.GetSampleID();
    int selection_ID = event.GetSelectionID();
    double true_beam_interactingEnergy = event.GetTrueInteractingEnergy();
    double reco_beam_interactingEnergy = event.GetRecoInteractingEnergy();
    double true_beam_endP = event.GetTrueEndP();
    //const std::vector<double> & true_beam_daughter_startP
    //    = event.GetTrueDaughterStartPs();
    const std::vector<double> & reco_beam_incidentEnergies
        = event.GetRecoIncidentEnergies();
    double reco_beam_endZ = event.GetRecoEndZ();
    //const std::vector<int> & true_beam_daughter_PDG
    //    = event.GetTrueDaughterPDGs();
    const std::vector<double> & true_beam_traj_Z = event.GetTrueTrajZ();
    const std::vector<double> & true_beam_traj_KE = event.GetTrueTrajKE();
    const std::vector<int>    & true_beam_slices = event.GetTrueSlices();
    const std::vector<double> & true_beam_incidentEnergies
        = event.GetTrueIncidentEnergies();

    if (samples.find(sample_ID) == samples.end())
      continue;

    double end_energy = true_beam_interactingEnergy;
    if (fSliceMethod == "Traj") {
      end_energy = sqrt(true_beam_endP*true_beam_endP*1.e6 + 139.57*139.57) - 139.57;
    }
    else if (fSliceMethod == "E") {
      end_energy = sqrt(true_beam_endP*true_beam_endP*1.e6 + 139.57*139.57) - 139.57;
    }
    else if (fSliceMethod == "Alt") {
      int bin = fEndSlices->GetXaxis()->FindBin(true_beam_traj_Z.back());
      if (bin > 0) {
        end_energy = fMeans.at(bin);
      }
    }

    double beam_P = event.GetBeamInstP();
    int scale_bin = GetBeamBin(bin_edges, beam_P);
    //std::cout << beam_P*1.e3 << " " << scale_bin << std::endl;
    double scale = variations[scale_bin];

    bool is_signal = signal_sample_checks.at(sample_ID);
    ThinSliceSample * this_sample = 0x0;
    if (is_signal) {
      std::vector<ThinSliceSample> & samples_vec = samples[sample_ID][0];
      //Get the samples vec from the first beam energy bin
      bool found = false;
      for (size_t j = 1; j < samples_vec.size()-1; ++j) {
        ThinSliceSample & sample = samples_vec.at(j);
        if (sample.CheckInSignalRange(end_energy)) {

          all_scales.push_back(scale);
          found = true;
          sample_scales[sample_ID][j] += scale;
          nominal_samples[sample_ID][j] += 1.;
          this_sample = &sample;
          break;
        }
      }
      if (!found) {
        if (end_energy < samples_vec[1].RangeLowEnd()) {

          all_scales.push_back(scale);
          sample_scales[sample_ID][0] += scale;
          nominal_samples[sample_ID][0] += 1.;
          this_sample = &samples_vec[0];
        }
        else {
          all_scales.push_back(scale);
          this_sample = &samples_vec.back();
          sample_scales[sample_ID].back() += scale;
          nominal_samples[sample_ID].back() += 1.;
        }
      }
    }
    else {
      this_sample = &samples[sample_ID][0][0];
      sample_scales[sample_ID][0] += scale;
      nominal_samples[sample_ID][0] += 1.;
    }

    new_flux += scale; //1 or scaled
    double val = 0.;
    if (selection_ID == 4) {
      if (selected_hists[selection_ID]->FindBin(reco_beam_endZ) == 0) {
        val = selected_hists[selection_ID]->GetBinCenter(1);
      }
      else if (selected_hists[selection_ID]->FindBin(reco_beam_endZ) >
               selected_hists[selection_ID]->GetNbinsX()) {
        val = selected_hists[selection_ID]->GetBinCenter(
            selected_hists[selection_ID]->GetNbinsX());
      }
      else {
        val = reco_beam_endZ;
      }
    }
    else if (selection_ID > 4) {
      val = .5;
    }
    else if (reco_beam_incidentEnergies.size()) {
      for (size_t j = 0; j < reco_beam_incidentEnergies.size(); ++j) {
        incident_hist.Fill(reco_beam_incidentEnergies[j]);
      }
      if (selected_hists.find(selection_ID) != selected_hists.end()) {
        //if (selection_ID != 4 && selection_ID != 5 && selection_ID != 6) {
        if (selection_ID < 4) {
          double energy = reco_beam_interactingEnergy;
          if (fDoEnergyFix) {
            for (size_t k = 1; k < reco_beam_incidentEnergies.size(); ++k) {
              double deltaE = (reco_beam_incidentEnergies[k-1] -
                               reco_beam_incidentEnergies[k]);
              if (deltaE > fEnergyFix) {
                energy += deltaE; 
              }
            }
          }
          if (selected_hists[selection_ID]->FindBin(energy) == 0) {
            val = selected_hists[selection_ID]->GetBinCenter(1);
          }
          else if (selected_hists[selection_ID]->FindBin(energy) >
                   selected_hists[selection_ID]->GetNbinsX()) {
            val = selected_hists[selection_ID]->GetBinCenter(
                selected_hists[selection_ID]->GetNbinsX());
          }
          else {
            val = energy;
          }
        }
      }
    }
    else {
      val = selected_hists[selection_ID]->GetBinCenter(1);
    }

    selected_hists[selection_ID]->Fill(val, scale);
    this_sample->AddVariedFlux(scale);
    std::vector<double> good_true_incEnergies = MakeTrueIncidentEnergies(
        true_beam_traj_Z, true_beam_traj_KE, true_beam_slices,
        true_beam_incidentEnergies);
    this_sample->AddIncidentEnergies(good_true_incEnergies, scale);
    if (norm_to_data_beam_P) {
      beam_fluxes[GetBeamBin(beam_energy_bins, beam_P)] += scale;
    }
  }

  double mean = 0.;
  for (auto s : all_scales) mean += s;
  mean /= all_scales.size();
  std::cout << "Mean scale: " << mean << std::endl;

  if (!norm_to_data_beam_P) {
    for (auto it = sample_scales.begin(); it != sample_scales.end(); ++it) {
      for (size_t i = 0; i < it->second.size(); ++i) {
        if (it->second[i] > 0.) {
          it->second[i] /= nominal_samples[it->first][i];
        }
        else {
          it->second[i] = 1.;
        }
        it->second[i] *= (flux/new_flux);
      }
    }

    incident_hist.Scale(flux/new_flux);
    for (auto it = selected_hists.begin(); it != selected_hists.end(); ++it) {
      it->second->Scale(flux/new_flux);
    }
  }
}

std::pair<double, size_t> protoana::AbsCexDriver::CalculateChi2(
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    ThinSliceDataSet & data_set) {

  double chi2 = 0.;
  double alt_chi2 = 0.;
  size_t nPoints = 0;
  size_t alt_nPoints = 0;

  std::map<int, TH1 *> & selected_data_hists = data_set.GetSelectionHists();
  // double total_data = 0., total_mc = 0.; // unused
  // double data_integral = 0.; // unused
  for (auto it = selected_data_hists.begin();
       it != selected_data_hists.end(); ++it) {
    TH1D * data_hist = (TH1D*)it->second;
    int selection_ID = it->first;
    if (data_hist->GetBinContent(0) > 0.) {
      std::cout << "Warning: underflow bin of " << selection_ID <<
                   " has " << data_hist->GetBinContent(0) << " events" <<
                   std::endl;
    }
    else if (data_hist->GetBinContent(data_hist->GetNbinsX()+1) > 0.) {
      std::cout << "Warning: overflow bin of " << selection_ID <<
                   " has " << data_hist->GetBinContent(data_hist->GetNbinsX()+1) <<
                   " events" <<
                   std::endl;
    }


    // data_integral += data_hist->Integral(); // unused
    int start = (fSkipFirstLast ? 2 : 1);
    int end = data_hist->GetNbinsX();
    if (fSkipFirstLast) --end;
    for (int i = start; i <= end; ++i) {
      double data_val = data_hist->GetBinContent(i);
      //double data_err = data_hist->GetBinError(i);

      double mc_val = 0.;
      double mc_sumw2 = 0.;
      //Go through all the samples and get the values from mc
      for (auto it2 = samples.begin(); it2 != samples.end(); ++it2) {
        std::vector<std::vector<ThinSliceSample>> & samples_vec_2D = it2->second;
        for (size_t j = 0; j < samples_vec_2D.size(); ++j) {
          std::vector<ThinSliceSample> & samples_vec = samples_vec_2D[j];
          for (size_t k = 0; k < samples_vec.size(); ++k) {
            ThinSliceSample & sample = samples_vec[k];
            mc_val += sample.GetSelectionHist(selection_ID)->GetBinContent(i);
            mc_sumw2 += std::pow(
                sample.GetSelectionHist(selection_ID)->GetBinError(i), 2);

            TH1D * mc_hist = (TH1D*)sample.GetSelectionHist(selection_ID);
            if (mc_hist->GetBinContent(0) > 0.) {
              std::cout << "Warning: underflow bin of " << selection_ID <<
                           " has " << mc_hist->GetBinContent(0) << " events" <<
                           std::endl;
            }
            else if (mc_hist->GetBinContent(mc_hist->GetNbinsX()+1) > 0.) {
              std::cout << "Warning: overflow bin of " << selection_ID <<
                           " has " << mc_hist->GetBinContent(mc_hist->GetNbinsX()+1) <<
                           " events" <<
                           std::endl;
            }
          }
        }
      }

      /*
      if (selection_ID != 6 && selection_ID != 5) {
        alt_chi2 += (std::pow((data_val - mc_val), 2) / mc_val);
        ++alt_nPoints;
      }*/


      double sigma_squared = (mc_val > 1.e-7 ? mc_sumw2/(mc_val*mc_val) : 1.);
      double beta = 1.;
      if (!fMultinomial) {
        beta = .5*((1. - mc_val*sigma_squared) +
                   sqrt(std::pow((1. - mc_val*sigma_squared), 2) +
                        4.*data_val*sigma_squared));
      }
      else {
        beta = .5 + .5*sqrt(1 + 4*data_val*sigma_squared);
      }

      /// Skip any bins with data == 0
      //
      //See PDG Stat Review:
      //https://pdg.lbl.gov/2018/reviews/rpp2018-rev-statistics.pdf
      //Page 6
      //
      //std::cout << mc_val << " " << data_val << " " << 2*data_val*std::log(data_val/mc_val) << std::endl;
      if (data_val > 1.e-7) {
        //chi2 += 2*data_val*std::log(data_val/mc_val)/* +
        //          (fBarlowBeeston ? ((beta - 1.)*(beta - 1.)/sigma_squared) : 0.)*/;
        chi2 += 2*data_val*std::log(data_val/(mc_val*(fBarlowBeeston ? beta : 1.))) +
                  (fBarlowBeeston ? ((beta - 1.)*(beta - 1.)/sigma_squared) : 0.);

        if (std::find(fToSkip.begin(), fToSkip.end(), selection_ID) ==
            fToSkip.end()) {
        //if (selection_ID != 6 && selection_ID != 5) {
          //std::cout << selection_ID << " " << mc_val << " " << beta << " " <<
          //             beta*mc_val << std::endl;
          double alt_chi2_cont = 2.*((fBarlowBeeston ? beta : 1.)*mc_val - data_val +
                      data_val*std::log(data_val/(mc_val*(fBarlowBeeston ? beta : 1.)))) +
                      (fBarlowBeeston ? ((beta - 1.)*(beta - 1.)/sigma_squared) : 0.);
          if (mc_val < 1.e-7) {
            std::cout << "Warning: " << selection_ID << " " << i << " mc_val "
                      << mc_val << std::endl;
            std::cout << "\t" << alt_chi2_cont << " "  << beta << " "
                      << data_val << " " << sigma_squared << std::endl;
          }
          alt_chi2 += 2.*((fBarlowBeeston ? beta : 1.)*mc_val - data_val +
                      data_val*std::log(data_val/(mc_val*(fBarlowBeeston ? beta : 1.)))) +
                      (fBarlowBeeston ? ((beta - 1.)*(beta - 1.)/sigma_squared) : 0.);
          //std::cout << "\t" << 2.*((fBarlowBeeston ? beta : 1.)*mc_val - data_val) << " " <<
          //             2.*data_val*std::log(data_val/(mc_val*(fBarlowBeeston ? beta : 1.)))
          //          << " " << (fBarlowBeeston ? ((beta - 1.)*(beta - 1.)/sigma_squared) : 0.)
          //          << std::endl;
          ++alt_nPoints;
        }
      }
      if (std::isnan(chi2) && fMultinomial) {
        std::cout << "Warning: " << selection_ID << " " << i << " data_val " <<
                     data_val << " mc_val " << mc_val << " isnan" << std::endl;
      }
      if (std::isnan(alt_chi2) && !fMultinomial) {
        std::cout << "Warning: " << selection_ID << " " << i << " data_val " <<
                     data_val << " mc_val " << mc_val << " chi2 isnan" << std::endl;
      }
      ++nPoints;
      // total_mc += mc_val; // unused
      // total_data += data_val; // unused
    }
    //std::cout << "Totals: " << total_data << " " << total_mc << std::endl;
  }
  //std::cout << std::endl;

  //if (total_data != total_mc) {
  //  std::cout << "WARNING: data does not match mc " << total_data << " " << total_mc << std::endl;
  //}

  //if (abs(chi2) < 1.e-7) chi2 = 0.;
  //std::cout << "totals (data, mc): " << total_data << " " << total_mc << std::endl;
  //if (chi2 < 0.) {
  //  std::cout << "Warning: chi2 < 0. " << std::setprecision(20) << total_data << " " << total_mc <<
  //               " " << chi2 << std::endl;
  //  std::cout << "Data Integral " << data_integral << std::endl;
  //}

  //-1 for multinomial
  std::pair<double, size_t> chi2_nPoints = {chi2, nPoints-1};
  std::pair<double, size_t> alt_chi2_nPoints = {alt_chi2, nPoints};
  return (fMultinomial ?  chi2_nPoints : alt_chi2_nPoints);
}

void protoana::AbsCexDriver::CompareSelections(
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    ThinSliceDataSet & data_set, TFile & output_file,
    std::vector<std::pair<int, int>> plot_style,
    bool plot_rebinned,
    bool post_fit, int nPars,
    TDirectory * plot_dir) {

  plot_dir->cd();
  //output_file.cd();
  std::map<int, TH1*> data_hists
      = (plot_rebinned ?
         data_set.GetRebinnedSelectionHists() :
         data_set.GetSelectionHists());
  for (auto it = data_hists.begin(); it != data_hists.end(); ++it) {
    int selection_ID = it->first;
    TH1D * data_hist = (TH1D*)it->second;
    data_hist->SetLineColor(kBlack);
    data_hist->SetMarkerColor(kBlack);
    data_hist->SetMarkerStyle(20);

    std::string canvas_name_no_data = "c";
    canvas_name_no_data += (post_fit ? "PostFit" : "Nominal") +
                   data_set.GetSelectionName(selection_ID) +
                   "NoData";
    TCanvas cSelectionNoData(canvas_name_no_data.c_str(), "");
    cSelectionNoData.SetTicks();

    std::string canvas_name = "c";
    canvas_name += (post_fit ? "PostFit" : "Nominal") +
                   data_set.GetSelectionName(selection_ID);
    TCanvas cSelection(canvas_name.c_str(), "");
    cSelection.SetTicks();

    std::map<int, std::vector<TH1D *>> temp_hists;
    for (auto it2 = samples.begin(); it2 != samples.end(); ++it2) {
      temp_hists[it2->first] = std::vector<TH1D *>();
      std::vector<ThinSliceSample> & vec = it2->second[0];
      for (size_t i = 0; i < vec.size(); ++i) {
        temp_hists[it2->first].push_back((TH1D*)(
            plot_rebinned ?
            vec[i].GetRebinnedSelectionHist(selection_ID) :
            vec[i].GetSelectionHist(selection_ID))->Clone());
      }
      for (size_t i = 1; i < it2->second.size(); ++i) {
        for (size_t j = 0; j < it2->second[i].size(); ++j) {
          temp_hists[it2->first][j]->Add((TH1D*)(
              plot_rebinned ?
              it2->second[i][j].GetRebinnedSelectionHist(selection_ID) :
              it2->second[i][j].GetSelectionHist(selection_ID)));
          }
      }
    }
    //Build the mc stack
    std::string stack_name = (post_fit ? "PostFit" : "Nominal") +
                             data_set.GetSelectionName(selection_ID) +
                             "Stack";
    THStack mc_stack(stack_name.c_str(), "");
    TLegend leg;
    std::vector<TH1D*> temp_vec;
    size_t iColor = 0;
    //need to add second loop with temp hists
    for (auto it2 = temp_hists.begin(); it2 != temp_hists.end(); ++it2) {
      for (size_t i = 0; i < it2->second.size(); ++i) {
        TH1D * sel_hist = it2->second.at(i);
        std::pair<int, int> color_fill = GetColorAndStyle(iColor, plot_style);
        sel_hist->SetFillColor(color_fill.first);
        sel_hist->SetFillStyle(color_fill.second);
        sel_hist->SetLineColor(kBlack);
        mc_stack.Add(sel_hist);
        temp_vec.push_back(sel_hist);
        ++iColor;
      }
    }
    mc_stack.Write();

    for (auto it = temp_vec.rbegin(); it != temp_vec.rend(); ++it) {
      leg.AddEntry(*it);
    }

    cSelectionNoData.cd();
    mc_stack.Draw("hist");
    std::string title_no_data = ";";
    title_no_data += data_hist->GetXaxis()->GetTitle();
    mc_stack.SetTitle(title_no_data.c_str());
    mc_stack.GetHistogram()->SetTitleSize(.04, "X");
    mc_stack.Draw("hist");
    if (it == data_hists.begin())
      leg.Write("leg_no_data");
    cSelectionNoData.Write();
    leg.Draw("same");

    leg.AddEntry(data_hist, "Data");
    
    std::pair<double, size_t> chi2 = CalculateChi2(samples, data_set);
    if (chi2.first < 0. && chi2.first > -1.e7) {
      chi2.first = 0.;
    }
    TString chi2_str;
    chi2_str.Form("#chi^{2} = %.2f", chi2.first);
    leg.AddEntry((TObject*)0x0, chi2_str, "");

    cSelection.cd();
    std::string title = data_set.GetSelectionName(selection_ID);
    title += ";";
    title += data_hist->GetXaxis()->GetTitle();
    mc_stack.SetTitle(title.c_str());

    mc_stack.Draw("hist");
    double max_mc = mc_stack.GetHistogram()->GetMaximum();
    int max_data_bin = data_hist->GetMaximumBin();
    double max_data = data_hist->GetBinContent(max_data_bin) +
                      data_hist->GetBinError(max_data_bin);
    mc_stack.SetMaximum((max_data > max_mc ? max_data : max_mc));
    mc_stack.Draw("hist");
    data_hist->Draw("e1 same");
    leg.Draw("same");

    cSelection.Write();

    //Get the full incident hist from stack
    TList * l = (TList*)mc_stack.GetHists();
    TH1D * hMC = (TH1D*)l->At(0)->Clone();
    for (int i = 1; i < l->GetSize(); ++i) {
      hMC->Add((TH1D*)l->At(i));
    }

    std::string ratio_name = data_set.GetSelectionName(selection_ID) + "Ratio" +
                             (post_fit ? "PostFit" : "Nominal");
    TH1D * hRatio
        = (TH1D*)data_hist->Clone(ratio_name.c_str());
    hRatio->Divide(hMC);
    hRatio->Write(); 
    std::string total_name = (post_fit ? "PostFit" : "Nominal") +
                             data_set.GetSelectionName(selection_ID) +
                             "Total";
    hMC->Write(total_name.c_str());

    canvas_name += "Ratio";
    TCanvas cRatio(canvas_name.c_str(), "");
    cRatio.SetTicks();
    TPad p1((canvas_name + "pad1").c_str(), "", 0, 0.3, 1., 1.);
    p1.SetBottomMargin(0.1);
    p1.Draw();
    p1.cd();
    mc_stack.Draw("hist");
    mc_stack.GetHistogram()->SetTitle("Abs;;");
    for (int i = 1; i < mc_stack.GetHistogram()->GetNbinsX(); ++i) {
      hRatio->GetXaxis()->SetBinLabel(
          i, mc_stack.GetHistogram()->GetXaxis()->GetBinLabel(i));
      mc_stack.GetHistogram()->GetXaxis()->SetBinLabel(i, "");
    }
    mc_stack.Draw("hist");
    data_hist->Draw("e1 same");

    cRatio.cd();
    TPad p2((canvas_name + "pad2").c_str(), "", 0, 0, 1., 0.3);
    p2.SetTopMargin(0.1);
    p2.SetBottomMargin(.2);
    p2.Draw();
    p2.cd();
    hRatio->Sumw2();
    std::string r_title = "";
    hRatio->GetYaxis()->SetTitle("Data/MC");
    hRatio->SetTitleSize(.04, "XY");
    hRatio->Draw("ep");

    cRatio.Write();


    //Differences
    //Get the full incident hist from stack
    TH1D * hMC2 = (TH1D*)l->At(0)->Clone();
    for (int i = 1; i < l->GetSize(); ++i) {
      hMC2->Add((TH1D*)l->At(i));
    }

    std::string diff_name = data_set.GetSelectionName(selection_ID) + "Diff" +
                             (post_fit ? "PostFit" : "Nominal");
    TH1D * hDiff
        = (TH1D*)data_hist->Clone(diff_name.c_str());
    hMC2->Scale(-1.);
    hDiff->Add(hMC2);
    hMC2->Scale(-1.);
    hDiff->Divide(hMC2);
    hDiff->Write(); 

    canvas_name += "Diff";
    TCanvas cDiff(canvas_name.c_str(), "");
    cDiff.SetTicks();
    cDiff.cd();
    hDiff->GetYaxis()->SetTitle("r");
    hDiff->SetTitleSize(.04, "XY");
    hDiff->Draw("ep");

    cDiff.Write();

    ///Write in NoStacks here

    THStack full_mc_stack((stack_name + "Full").c_str(), "");
    //TLegend leg;
    size_t iColorFull = 0;
    //need to add second loop with temp hists
    for (auto it2 = temp_hists.begin(); it2 != temp_hists.end(); ++it2) {
      TH1D * sel_hist = it2->second.at(0);
      std::pair<int, int> color_fill = GetColorAndStyle(iColorFull, plot_style);
      sel_hist->SetFillColor(color_fill.first);
      sel_hist->SetFillStyle(color_fill.second);
      sel_hist->SetLineColor(kBlack);


      for (size_t i = 1; i < it2->second.size(); ++i) {
        sel_hist->Add(it2->second.at(i));
      }
      temp_vec.push_back(sel_hist);
      full_mc_stack.Add(sel_hist);
      ++iColorFull;
    }
    full_mc_stack.Write();
  }

  /* total and total_muons set but not used
  double total_muons = 0.; // unused
  double total = 0.;
  for (auto it = samples.begin(); it != samples.end(); ++it) {
    for (size_t i = 0; i < it->second.size(); ++i) {
      for (size_t j = 0; j < it->second[i].size(); ++j) {
        if (it->first == 5) {
          total_muons += it->second[i][j].GetVariedFlux();
        }
        total += it->second[i][j].GetVariedFlux();
      }
    }
  }
  */
}

void protoana::AbsCexDriver::GetCurrentHists(
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    ThinSliceDataSet & data_set,
    std::map<int, std::vector<TH1*>> & throw_hists,
    bool plot_rebinned) {
  
  //Use data set as a template
  std::map<int, TH1*> data_hists
      = (plot_rebinned ?
         data_set.GetRebinnedSelectionHists() :
         data_set.GetSelectionHists());

  for (auto it = data_hists.begin(); it != data_hists.end(); ++it) {
    std::vector<double> bins;
    std::string name = data_set.GetSelectionName(it->first) + "Throw" +
                       std::to_string(throw_hists[it->first].size());

    TH1D * temp_hist = (TH1D*)it->second->Clone(name.c_str());
    temp_hist->Reset();
    throw_hists[it->first].push_back(temp_hist);
  }

  //Iterate over samples
  for (auto it = throw_hists.begin(); it != throw_hists.end(); ++it) {
    int selection_ID = it->first;
    for (auto it2 = samples.begin(); it2 != samples.end(); ++it2) {
      std::vector<std::vector<ThinSliceSample>> & samples_vec_2D = it2->second;
      for (size_t j = 0; j < samples_vec_2D.size(); ++j) {
        std::vector<ThinSliceSample> & samples_vec = samples_vec_2D[j];
        for (size_t k = 0; k < samples_vec.size(); ++k) {
          ThinSliceSample & sample = samples_vec[k];
          for (int i = 1; i <= it->second.back()->GetNbinsX(); ++i) {
            it->second.back()->AddBinContent(i,
                sample.GetSelectionHist(selection_ID)->GetBinContent(i));
          }
        }
      }
    }
  }
}

void protoana::AbsCexDriver::GetCurrentTruthHists(
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    std::map<int, std::vector<TH1*>> & throw_hists,
    std::map<int, std::vector<TH1*>> & throw_inc_hists,
    std::map<int, std::vector<TH1*>> & throw_xsec_hists,
    const std::vector<int> & incident_samples,
    const std::map<int, std::vector<double>> & signal_bins) {
  //Loop over the samples
  for (auto it = samples.begin(); it != samples.end(); ++it) {
    //Get the number of bins from the first entry of the beam energy bins
    std::vector<std::vector<ThinSliceSample>> & samples_vec_2D = it->second;
    size_t nBins = samples_vec_2D[0].size();
    std::string name = it->second[0][0].GetName() + "Throw" +
                       std::to_string(throw_hists[it->first].size());
    TH1D * temp_hist = new TH1D(name.c_str(), "", nBins, 0, nBins); 
    for (size_t i = 0; i < samples_vec_2D.size(); ++i) {
      std::vector<ThinSliceSample> & samples_vec = samples_vec_2D[i];
      for (size_t j = 0; j < it->second[i].size(); ++j) {
        temp_hist->AddBinContent(j+1, samples_vec[j].GetVariedFlux());
      }
    }
    throw_hists[it->first].push_back(temp_hist);
  }

  for (auto it = throw_inc_hists.begin(); it != throw_inc_hists.end(); ++it) {
    int s = it->first;
    auto & samples_vec_2D = samples[s];
    const std::vector<double> & bins = signal_bins.at(s);
    std::string name = samples_vec_2D[0][0].GetName();
    name += "IncidentThrow" +
             std::to_string(throw_inc_hists[it->first].size());
    TH1D * temp_inc_hist = new TH1D(name.c_str(), "", bins.size() - 1, &bins[0]); 
    

    name = samples_vec_2D[0][0].GetName();
    name += "XSecThrow" +
             std::to_string(throw_inc_hists[it->first].size());
    TH1D * temp_xsec_hist = new TH1D(name.c_str(), "", bins.size() - 1,
                                     &bins[0]);
    for (auto i_s : incident_samples) {
      auto & incident_vec_2D = samples[i_s];
      for (size_t i = 0; i < incident_vec_2D.size(); ++i) {
        for (size_t j = 0; j < incident_vec_2D[i].size(); ++j) {
          //if (fExtraOptions.get<std::string>("SliceMethod") == "E") {
          if (/*fExtraOptions.get<std::string>("SliceMethod")*/fSliceMethod == "E") {
            incident_vec_2D[i][j].FillESliceHist(*temp_inc_hist);
          }
          else {
            incident_vec_2D[i][j].FillHistFromIncidentEnergies(*temp_inc_hist);
          }
        }
      }
    }
    throw_inc_hists[s].push_back(temp_inc_hist);

    for (int i = 1; i <= temp_xsec_hist->GetNbinsX(); ++i) {
      temp_xsec_hist->SetBinContent(
          i, throw_hists[s].back()->GetBinContent(i+1));
    }
    temp_xsec_hist->Divide(temp_inc_hist);
    throw_xsec_hists[s].push_back(temp_xsec_hist);
  }
}

void protoana::AbsCexDriver::PlotThrows(
    ThinSliceDataSet & data_set, std::map<int, std::vector<TH1*>> & throw_hists,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    size_t nThrows,
    std::map<int, std::vector<TH1*>> & truth_throw_hists,
    std::map<int, std::vector<TH1*>> & truth_inc_hists,
    std::map<int, std::vector<TH1*>> & truth_xsec_hists,
    std::map<int, TH1*> & best_fit_incs,
    std::map<int, TH1*> & best_fit_xsecs,
    std::map<int, TH1*> & nominal_incs,
    std::map<int, TH1*> & nominal_xsecs,
    TFile & output_file, bool plot_rebinned,
    std::map<int, std::vector<double>> * sample_scales) {
  std::map<int, TH1*> data_hists
      = (plot_rebinned ?
         data_set.GetRebinnedSelectionHists() :
         data_set.GetSelectionHists());

  //Build best fit hists and get bins for covariance 
  std::map<int, TH1D*> best_fit_selection_hists;
  int nBins = 0;
  for (auto it = data_hists.begin(); it != data_hists.end(); ++it ) {
    TH1D * best_fit_hist = (TH1D*)it->second->Clone();
    best_fit_hist->Reset();
    for (auto it2 = samples.begin(); it2 != samples.end(); ++it2) {
      for (size_t i = 0; i < it2->second.size(); ++i) {
        for (size_t j = 0; j < it2->second[i].size(); ++j) {
          it2->second[i][j].SetFactorToBestFit();
          best_fit_hist->Add(
              (TH1D*)(plot_rebinned ?
                      it2->second[i][j].GetRebinnedSelectionHist(it->first) :
                      it2->second[i][j].GetSelectionHist(it->first)));
        }
      }
    }
    best_fit_selection_hists[it->first] = best_fit_hist;
    nBins += best_fit_hist->GetNbinsX();
  }

  TH2D selection_cov("SelectionCov", "", nBins, 0, nBins, nBins, 0, nBins);

  nBins = 0;
  std::map<int, size_t> sample_bins;
  for (auto it = samples.begin(); it != samples.end(); ++it) {
    nBins += it->second[0].size();
    sample_bins[it->first] = it->second[0].size();
  }

  std::map<int, std::vector<double>> best_fit_truth;
  std::map<int, std::vector<double>> best_fit_errs;

  for (auto it = samples.begin(); it != samples.end(); ++it) {
    best_fit_truth[it->first]
        = std::vector<double>(sample_bins[it->first], 0.);
    best_fit_errs[it->first]
        = std::vector<double>(sample_bins[it->first], 0.);
   
    for (size_t i = 0; i < sample_bins[it->first]; ++i) {
      double best_fit_val_i = 0.;
      for (size_t j = 0; j < it->second.size(); ++j) {
        best_fit_val_i += it->second[j][i].GetVariedFlux();
      }

      best_fit_truth[it->first][i] = best_fit_val_i;
    }
  }

  TH2D interaction_cov("interaction_cov", "", nBins, 0, nBins, nBins, 0, nBins);
  std::map<int, std::vector<double>> best_fit_inc_truth;
  std::map<int, std::vector<double>> best_fit_xsec_truth;
  std::map<int, std::vector<double>> best_fit_inc_errs;
  std::map<int, std::vector<double>> best_fit_xsec_errs;

  nBins = 0;
  std::map<int, size_t> xsec_bins;
  for (auto it = best_fit_incs.begin(); it != best_fit_incs.end(); ++it) {
    int s = it->first;
    nBins += it->second->GetNbinsX();
    xsec_bins[s] = it->second->GetNbinsX();

    best_fit_inc_truth[s] = std::vector<double>(xsec_bins[s], 0.);
    best_fit_xsec_truth[s] = std::vector<double>(xsec_bins[s], 0.);
    best_fit_inc_errs[s] = std::vector<double>(xsec_bins[s], 0.);
    best_fit_xsec_errs[s] = std::vector<double>(xsec_bins[s], 0.);
    
    for (size_t i = 0; i < xsec_bins[s]; ++i) {
      best_fit_inc_truth[s][i] = it->second->GetBinContent(i+1);
      best_fit_xsec_truth[s][i] = best_fit_xsecs[s]->GetBinContent(i+1);
    }
  }

  //TH2D incident_cov("incident_cov", "", nBins, 0, nBins, nBins, 0, nBins);
  TH2D xsec_cov("xsec_cov", "", nBins, 0, nBins, nBins, 0, nBins);

  for (size_t z = 0; z < nThrows; ++z) {
    int bin_i = 1;
    for (auto it = best_fit_selection_hists.begin();
         it != best_fit_selection_hists.end(); ++it) {
      TH1D * best_fit = it->second;
      int selection_ID = it->first;
      std::vector<TH1*> & temp_throws = throw_hists[selection_ID];
      for (int i = 1; i <= best_fit->GetNbinsX(); ++i) {
        double best_fit_val_i = best_fit->GetBinContent(i);
        int bin_j = 1;
        for (auto it2 = best_fit_selection_hists.begin();
             it2 != best_fit_selection_hists.end(); ++it2) {

          TH1D * best_fit_2 = it2->second;
          int selection_ID_2 = it2->first;
          std::vector<TH1*> & temp_throws_2 = throw_hists[selection_ID_2];
          for (int j = 1; j <= best_fit_2->GetNbinsX(); ++j) {
            double best_fit_val_j = best_fit_2->GetBinContent(j);
            double val = (best_fit_val_i - temp_throws[z]->GetBinContent(i))*
                         (best_fit_val_j - temp_throws_2[z]->GetBinContent(j));
            selection_cov.SetBinContent(
                bin_i, bin_j, (val/temp_throws.size() +
                               selection_cov.GetBinContent(bin_i, bin_j)));
            ++bin_j;
          }
        }
        ++bin_i;
      }
    }

    bin_i = 1;
    for (auto it = samples.begin(); it != samples.end(); ++it) {
      std::vector<TH1 *> throw_hists_i = truth_throw_hists[it->first];
     
      for (size_t i = 0; i < sample_bins[it->first]; ++i) {
        double best_fit_val_i = best_fit_truth[it->first][i];

        int bin_j = 1;
        for (auto it2 = samples.begin(); it2 != samples.end(); ++it2) {
          std::vector<TH1 *> throw_hists_j = truth_throw_hists[it2->first];
          for (size_t j = 0; j < sample_bins[it2->first]; ++j) {
            double best_fit_val_j = best_fit_truth[it2->first][j];

            double val
                = (throw_hists_i[z]->GetBinContent(i+1) - best_fit_val_i)*
                  (throw_hists_j[z]->GetBinContent(j+1) - best_fit_val_j);
            interaction_cov.SetBinContent(
                bin_i, bin_j,
                (interaction_cov.GetBinContent(bin_i, bin_j) +
                 val/throw_hists_i.size()));
            if (bin_i == bin_j && (z == nThrows - 1)) {
              best_fit_errs[it->first][i]
                  = sqrt(interaction_cov.GetBinContent(bin_i, bin_j));
            }
            ++bin_j;
          }
        }

        ++bin_i;
      }
    }

    bin_i = 1;
    for (auto it = truth_inc_hists.begin(); it != truth_inc_hists.end(); ++it) {
      //std::vector<TH1 *> inc_hists_i = it->second;
      std::vector<TH1 *> xsec_hists_i = truth_xsec_hists[it->first];

      for (size_t i = 0; i < xsec_bins[it->first]; ++i) {
        //double best_fit_inc_i = best_fit_inc_truth[it->first][i];
        double best_fit_xsec_i = best_fit_xsec_truth[it->first][i];

        int bin_j = 1;
        for (auto it2 = truth_inc_hists.begin(); it2 != truth_inc_hists.end();
             ++it2) {
          std::vector<TH1 *> xsec_hists_j = truth_xsec_hists[it2->first];
          for (size_t j = 0; j < xsec_bins[it2->first]; ++j) {
            double best_fit_xsec_j = best_fit_xsec_truth[it2->first][j];

            double val
                = (xsec_hists_i[z]->GetBinContent(i+1) - best_fit_xsec_i)*
                  (xsec_hists_j[z]->GetBinContent(j+1) - best_fit_xsec_j);
            xsec_cov.SetBinContent(
                bin_i, bin_j,
                (xsec_cov.GetBinContent(bin_i, bin_j) +
                 val/nThrows));
            if (bin_i == bin_j && (z == nThrows - 1)) {
              best_fit_xsec_errs[it->first][i]
                  = sqrt(xsec_cov.GetBinContent(bin_i, bin_j));
            }
          }
        }
      }
    }
  }


  output_file.cd("Throws");
  selection_cov.Write();
  interaction_cov.Write();
  xsec_cov.Write();

  int bin_count = 0;
  for (auto it = data_hists.begin(); it != data_hists.end(); ++it) {
    int selection_ID = it->first;
    std::vector<TH1*> hists = throw_hists.at(selection_ID);

    std::string canvas_name = "cThrow" +
                              data_set.GetSelectionName(selection_ID);
    TCanvas cThrow(canvas_name.c_str(), "");
    cThrow.SetTicks();

    std::string name = "Throw" + data_set.GetSelectionName(selection_ID);
    auto data_hist = it->second;
    std::vector<double> xs, xs_width;
    std::vector<double> ys, errs;
    for (int i = 1;
         i <= best_fit_selection_hists[it->first]->GetNbinsX(); ++i) {
      ys.push_back(
          best_fit_selection_hists[it->first]->GetBinContent(i));
      errs.push_back(
          sqrt(selection_cov.GetBinContent(bin_count+i, bin_count+i)));
      xs.push_back(data_hist->GetBinCenter(i));
      xs_width.push_back(data_hist->GetBinWidth(i)/2.);
    } 

    TGraphAsymmErrors throw_gr(data_hist->GetNbinsX(),
                               &xs[0], &ys[0], 
                               &xs_width[0], &xs_width[0], &errs[0], &errs[0]);

    throw_gr.SetFillStyle(3144);
    throw_gr.SetFillColor(kRed);
    throw_gr.Draw("a2");
    data_hist->Draw("same e1");
    output_file.cd("Throws");
    cThrow.Write();

    bin_count += data_hist->GetNbinsX();
  }

  bin_count = 0;
  for (auto it = truth_throw_hists.begin(); it != truth_throw_hists.end(); ++it) {
    int sample_ID = it->first;

    std::vector<double> xs, xs_width;
    for (size_t i = 0; i < sample_bins[it->first]; ++i) {
      xs.push_back(i + 0.5);
      xs_width.push_back(.5);
    }

    std::string name = "hNominal" + samples[sample_ID][0][0].GetName();
    TH1D temp_nominal(name.c_str(), "", xs.size(), 0, xs.size());
    std::vector<std::vector<ThinSliceSample>> & samples_vec_2D
        = samples[sample_ID];
    for (size_t i = 0; i < samples_vec_2D.size(); ++i) {
      for (size_t j = 0; j < samples_vec_2D[i].size(); ++j) {
        temp_nominal.AddBinContent(j+1, samples_vec_2D[i][j].GetNominalFlux());
      }
    }

    double max = -999.;
    for (size_t i = 0; i < sample_bins[it->first]; ++i) {
      if ((best_fit_truth[sample_ID][i] + best_fit_errs[sample_ID][i]) > max)
        max = (best_fit_truth[sample_ID][i] + best_fit_errs[sample_ID][i]);

      if (temp_nominal.GetBinContent(i+1) > max)
        max = temp_nominal.GetBinContent(i+1);
    }

    output_file.cd("Throws");
    std::string canvas_name = "cTruthThrow" + samples[sample_ID][0][0].GetName();
    TCanvas cThrow(canvas_name.c_str(), "");
    cThrow.SetTicks();
    TGraphAsymmErrors throw_gr(xs.size(),
                                &xs[0], &best_fit_truth[it->first][0], 
                                &xs_width[0], &xs_width[0],
                                &best_fit_errs[it->first][0],
                                &best_fit_errs[it->first][0]);
    throw_gr.SetFillStyle(3144);
    throw_gr.SetFillColor(kRed);
    throw_gr.SetMinimum(0.);
    throw_gr.SetMaximum(1.5*max);
    throw_gr.Draw("a2");
    throw_gr.Draw("p");

    temp_nominal.SetMarkerColor(kBlue);
    temp_nominal.SetMarkerStyle(20);
    temp_nominal.Draw("same p");

    TLegend leg;
    leg.AddEntry(&throw_gr, "Throws", "lpf");
    leg.AddEntry(&temp_nominal, "Nominal", "p");

    if (sample_scales) {
      name = "hVaried" + samples[sample_ID][0][0].GetName();
      TH1D * temp_varied = (TH1D*)temp_nominal.Clone(name.c_str());
      for (size_t i = 0; i < xs.size(); ++i) {
        temp_varied->SetBinContent(
            i+1, temp_varied->GetBinContent(i+1)*(*sample_scales)[sample_ID][i]);
      }
      temp_varied->SetMarkerColor(kBlack);
      temp_varied->SetMarkerStyle(20);
      temp_varied->Draw("same p");
      leg.AddEntry(temp_varied, "Fake Data", "p");
    }

    leg.Draw();
    cThrow.Write();

    bin_count += xs.size();
  }

  for (auto it = best_fit_xsec_truth.begin(); it != best_fit_xsec_truth.end();
       ++it) {
    int sample_ID = it->first;

    std::vector<double> xs, xs_width;
    for (size_t i = 0; i < xsec_bins[sample_ID]; ++i) {
      xs.push_back(i + 0.5);
      xs_width.push_back(.5);
    }

    output_file.cd("Throws");
    std::string canvas_name = "cXSecThrow" + samples[sample_ID][0][0].GetName();
    TCanvas cThrow(canvas_name.c_str(), "");
    cThrow.SetTicks();
    TGraphAsymmErrors throw_gr(xs.size(),
                               &xs[0], &best_fit_xsec_truth[it->first][0], 
                               &xs_width[0], &xs_width[0],
                               &best_fit_xsec_errs[it->first][0],
                               &best_fit_xsec_errs[it->first][0]);
    throw_gr.SetFillStyle(3144);
    throw_gr.SetFillColor(kRed);
    throw_gr.SetMinimum(0.);
    throw_gr.Draw("a2");
    throw_gr.Draw("p");

    cThrow.Write();
  }
}

//Can Be removed
int protoana::AbsCexDriver::RecalculateSelectionID(
    const ThinSliceEvent & event,
    double C_cal,
    TProfile * prot_template) {

  if (event.GetSelectionID() > 3) {
    return event.GetSelectionID();
  }

  //Look for charged pions
  //bool has_pi0_shower = false;
  const std::vector<double> & track_scores = event.GetRecoDaughterTrackScores();
  for (size_t i = 0; i < track_scores.size(); ++i) {

    if (track_scores[i] > 0.3) {
      //Here: recalculate dEdX and truncated mean dEdX
      std::vector<double> new_dEdX, new_res_range;
      const std::vector<double> & calibrated_dQdX
          = event.GetRecoDaughterTrackdQdXs()[i];
      const std::vector<double> & daughter_EField
          = event.GetRecoDaughterEFields()[i];
      const std::vector<double> & res_range
          =event.GetRecoDaughterTrackResRanges()[i];
      for (size_t j = 0; j < calibrated_dQdX.size(); ++j) {
      //std::cout << C_cal << " " << (calibrated_dQdX)[j] << " " << fBetaP << " "
      //          << fRho << " " << daughter_EField[j] << " " << fWion << " "
      //          << fAlpha << std::endl;
        if (calibrated_dQdX[j] < 0)
          continue;
          //std::cout << "Warning dqdx < 0: " << calibrated_dQdX[j] << std::endl;
        double dedx = C_cal;//(1./*/(C_cal)*/);
        dedx *= (calibrated_dQdX)[j];
        dedx *= (fBetaP / ( fRho * (daughter_EField)[j] ) * fWion);
        dedx = exp(dedx);
        dedx -= fAlpha;
        dedx *= ((fRho*(daughter_EField)[j])/fBetaP);
        new_dEdX.push_back(dedx);
        new_res_range.push_back(res_range[j]);

        //std::cout << "Added " << dedx << std::endl;
      }

      double truncated_mean = TruncatedMean(new_dEdX);
      //std::cout << i << " trunc mean: " << truncated_mean << " "; 

      if (truncated_mean < 2.8 && truncated_mean > 0.5) {
       // std::cout << std::endl;
        return 3;
      }
      else if (truncated_mean > 2.8 && truncated_mean < 3.4) {
        std::pair<double, int> pid_chi2_ndof
            = fTrackUtil.Chi2PID(
                new_dEdX, new_res_range, prot_template);
       // std::cout <<  pid_chi2_ndof.first/pid_chi2_ndof.second;
        if (pid_chi2_ndof.second > 0 &&
            pid_chi2_ndof.first/pid_chi2_ndof.second > 70.) {
            //if (event.GetSelectionID() != 3) {
            //  std::cout << i << " chi2: " << pid_chi2_ndof.first/pid_chi2_ndof.second << " " << truncated_mean << std::endl;
            //  std::cout << "\t" << new_dEdX.size() << " " << new_res_range.size() << std::endl;
            //}
 //         std::cout << std::endl;
          return 3;
        }
      }
   //   std::cout << std::endl;
    }
  }

  if (event.GetHasPi0Shower()) return 2;

  return 1;
}

double protoana::AbsCexDriver::TruncatedMean(const std::vector<double> & dEdX) {
  if (!dEdX.size()) return -999.;
  std::vector<double> temp_dEdX;
  for (auto d : dEdX) {
    if (d < 0) continue;
    temp_dEdX.push_back(d);
  }
  std::sort(temp_dEdX.begin(), temp_dEdX.end());
  size_t low = std::rint(temp_dEdX.size()*0.16);
  size_t high = std::rint(temp_dEdX.size()*0.84);
  std::vector<double> temp;
  for (size_t i = low; i < high; ++i) {
    temp.push_back(temp_dEdX[i]);
  }
  if (temp.empty()) return -999.;
  return std::accumulate(temp.begin(), temp.end(), 0.0)/temp.size();
}

std::vector<double> protoana::AbsCexDriver::MakeTrueIncidentEnergies(
    const std::vector<double> & true_beam_traj_Z,
    const std::vector<double> & true_beam_traj_KE,
    const std::vector<int> & true_beam_slices,
    const std::vector<double> & true_beam_incidentEnergies) {

  std::vector<double> results;
  if (fSliceMethod == "Traj") {
    double next_slice_z = fTrajZStart;
    int next_slice_num = 0;
    for (size_t j = 1; j < true_beam_traj_Z.size() - 1; ++j) {
      double z = true_beam_traj_Z[j];
      double ke = true_beam_traj_KE[j];

      if (z < fTrajZStart) {
        continue;
      }

      if (z >= next_slice_z) {
        double temp_z = true_beam_traj_Z[j-1];
        double temp_e = true_beam_traj_KE[j-1];
        
        while (next_slice_z < z && next_slice_num < fSliceCut) {
          double sub_z = next_slice_z - temp_z;
          double delta_e = true_beam_traj_KE[j-1] - ke;
          double delta_z = z - true_beam_traj_Z[j-1];
          temp_e -= (sub_z/delta_z)*delta_e;
          results.push_back(temp_e);
          temp_z = next_slice_z;
          next_slice_z += fPitch;
          ++next_slice_num;
        }
      }
    }
  }
  else if (fSliceMethod == "Default") {
    for (size_t j = 0; j < true_beam_incidentEnergies.size(); ++j) {
      int slice = true_beam_slices[j]; 
      if (slice > fSliceCut) continue;
      results.push_back(true_beam_incidentEnergies[j]);
    }
  }
  else if (fSliceMethod == "Alt") {
    int bin = fEndSlices->GetXaxis()->FindBin(true_beam_traj_Z.back());
    for (int i = 1; i <= bin; ++i) {
      results.push_back(fMeans.at(i));
    }
  }

  return results;
}

int protoana::AbsCexDriver::GetBeamBin(
    const std::vector<double> & beam_energy_bins,
    const double & true_beam_startP,
    bool restrict_P) {
  int bin = -1;
  for (size_t j = 1; j < beam_energy_bins.size(); ++j) {
    if ((beam_energy_bins[j-1] <= 1.e3*true_beam_startP) &&
        (1.e3*true_beam_startP < beam_energy_bins[j])) {
      bin = j - 1;
      break;
    }
  }
  if (bin == -1 && !restrict_P) {
    std::string message = "Could not find beam energy bin for " +
                          std::to_string(true_beam_startP);
    throw std::runtime_error(message);
  }
  return bin;
}

void protoana::AbsCexDriver::ConstructCovariances(
    const std::vector<ThinSliceEvent> & events,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & nominal_samples,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & covariance_samples,
    const std::map<int, bool> & signal_sample_checks,
    std::vector<double> & beam_energy_bins,
    std::map<int, double> & nominal_fluxes,
    std::map<int, std::vector<std::vector<double>>> & fluxes_by_sample,
    const std::map<int, std::vector<double>> & signal_pars,
    const std::map<int, double> & flux_pars,
    const std::map<std::string, ThinSliceSystematic> & syst_pars,
    const std::map<std::string, ThinSliceSystematic> & g4rw_pars,
    bool fit_under_over, bool tie_under_over, bool use_beam_inst_P

    /*const std::vector<ThinSliceEvent> & events,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & nominal_samples,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & covariance_samples,
    const std::map<int, bool> & signal_sample_checks,
    std::map<int, double> & nominal_fluxes,
    std::map<int, std::vector<std::vector<double>>> & fluxes_by_sample,
    std::vector<double> & beam_energy_bins, bool use_beam_inst_P*/) {
  for (auto & routine : fCovarianceRoutines) {
    if (routine == "BeamShift") {
      CovarianceRoutineBeamShift(
          events, nominal_samples, covariance_samples, signal_sample_checks,
          beam_energy_bins, nominal_fluxes, fluxes_by_sample, signal_pars,
          flux_pars, syst_pars, g4rw_pars, fit_under_over, tie_under_over, use_beam_inst_P);
    }
  }
}

void protoana::AbsCexDriver::SetupBeamShiftCovRoutine(fhicl::ParameterSet & routine) {
  fBeamShiftCovP0 = routine.get<std::pair<double, double>>("P0");
  fBeamShiftCovP1 = routine.get<std::pair<double, double>>("P1");
  fBeamShiftCovP2 = routine.get<std::pair<double, double>>("P2");
  std::cout << fBeamShiftCovP0.first << " " << fBeamShiftCovP0.second << std::endl;
  std::cout << fBeamShiftCovP1.first << " " << fBeamShiftCovP1.second << std::endl;
  std::cout << fBeamShiftCovP2.first << " " << fBeamShiftCovP2.second << std::endl;

  fBeamShiftCovUseInput = routine.get<bool>("UseInput", false);
  if (fBeamShiftCovUseInput) {
    fBeamShiftInputFileName = routine.get<std::string>("InputName");
  }
    /*
    auto * beam_shift_file = TFile::Open(beam_shift_input_file.c_str());

    fBeamShiftInputCentrals = (TVectorD*)beam_shift_file->Get("centrals");
    //fBeamShiftInputCentrals->SetDirectory(0);

    fBeamShiftInputCov = (TMatrixD*)beam_shift_file->Get("covariance");
    //fBeamShiftInputCov->SetDirectory(0);

    fBeamShiftInputChol = TDecompChol(*fBeamShiftInputCov);
    fBeamShiftInputChol.Decompose();
    //fBeamShiftInputChol.SetDirectory(0);

    beam_shift_file->Close();
  }*/


  fBeamShiftCovOutput = routine.get<std::string>("OutputName");
  fNCovarianceGens = routine.get<size_t>("NCovarianceGens");

  fBeamShiftNominalP0 = routine.get<double>("NominalP0");
  fBeamShiftNominalP1 = routine.get<double>("NominalP1");
  fBeamShiftNominalP2 = routine.get<double>("NominalP2");
}

void protoana::AbsCexDriver::CovarianceRoutineBeamShift(
    const std::vector<ThinSliceEvent> & events,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & nominal_samples,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & new_samples,
    const std::map<int, bool> & signal_sample_checks,
    std::vector<double> & beam_energy_bins,
    std::map<int, double> & nominal_fluxes,
    std::map<int, std::vector<std::vector<double>>> & fluxes_by_sample,
    const std::map<int, std::vector<double>> & signal_pars,
    const std::map<int, double> & flux_pars,
    const std::map<std::string, ThinSliceSystematic> & syst_pars,
    const std::map<std::string, ThinSliceSystematic> & g4rw_pars,
    bool fit_under_over, bool tie_under_over, bool use_beam_inst_P) {

  fBeamShiftCovRoutineActive = true;
  if (fBeamShiftCovUseInput) {
    //std::string beam_shift_input_file = routine.get<std::string>("InputName");
    auto * beam_shift_file = TFile::Open(fBeamShiftInputFileName.c_str());

    fBeamShiftInputCentrals = (TVectorD*)beam_shift_file->Get("centrals");
    //fBeamShiftInputCentrals->SetDirectory(0);

    fBeamShiftInputCov = (TMatrixD*)beam_shift_file->Get("covariance");
    //fBeamShiftInputCov->SetDirectory(0);

    fBeamShiftInputChol = TDecompChol(*fBeamShiftInputCov);
    fBeamShiftInputChol.Decompose();
    fBeamShiftInputCovL = (TMatrixD*)fBeamShiftInputChol.GetU().Clone();
    fBeamShiftInputCovL->Transpose(fBeamShiftInputChol.GetU());
    //fBeamShiftInputChol.SetDirectory(0);

    beam_shift_file->Close();
  }

  //Nominal samples -- get total histograms
  std::vector<double> nominal_vals;
  bool push_back = true;
  for (auto it = nominal_samples.begin(); it != nominal_samples.end(); ++it) {
    for (size_t i = 0; i < it->second.size(); ++i) {
      for (size_t j = 0; j < it->second[i].size(); ++j) {
        const std::map<int, TH1*> & sel_hists
            = it->second[i][j].GetSelectionHists();
        int bin = 0;
        for (auto it2 = sel_hists.begin(); it2 != sel_hists.end(); ++it2) {
          if (it2->first == 6 || it2->first == 5) continue;
          //std::cout << it2->first << std::endl;

          for (int k = 1; k <= it2->second->GetNbinsX(); ++k) {
            if (push_back) {
              nominal_vals.push_back(it2->second->GetBinContent(k));
            }
            else {
              nominal_vals[bin] += it2->second->GetBinContent(k);
            }
            ++bin;
          }
        }
        if (push_back) push_back = false;
      }
    }
  }

  std::vector<double> mean_new_vals(nominal_vals.size());
  //std::vector<double> cov_vals(nominal_vals.size());

  //int a = 0;
  //for (auto & nv : nominal_vals) {
  //  std::cout << a << " " << nv << std::endl;
  //  ++a;
  //}

  TFile output_cov_file(fBeamShiftCovOutput.c_str(), "recreate");
  TH2D output_cov("cov_hist", "", nominal_vals.size(), 0, nominal_vals.size(),
                  nominal_vals.size(), 0, nominal_vals.size());
  TMatrixD output_cov_mat(nominal_vals.size(), nominal_vals.size());
  //for (size_t i = 0; i < cov_vals.size(); ++i) {
  //  for (size_t j = 0; j < cov_vals.size(); ++j) {
  //    output_cov_mat[i][j] = cov_vals[i]*cov_vals[j]/fNCovarianceGens;
  //    output_cov.SetBinContent(i+1, j+1,
  //                             cov_vals[i]*cov_vals[j]/fNCovarianceGens);
  //  }
  //}



  //
  //
  //Create TH2D/TMatrixD to hold relative differences
  //bin_i_j = (1/fNCovarianceGens)*sum((bin_i - bin_i_nom)*(bin_j - bin_j_nom)/(bin_i*bin_j))
  //int fNCovarianceGens = 2;
  std::cout << "Generating covariance for beam shift" << std::endl;
  std::vector<std::vector<double>> all_new_vals;
  for (size_t gen_i = 0; gen_i < fNCovarianceGens; ++gen_i) {

    /*
    fCurrentBeamShiftP0 = fBeamShiftCovP0.first + fBeamShiftCovP0.second*fRNG.Gaus();
    fCurrentBeamShiftP1 = fBeamShiftCovP1.first + fBeamShiftCovP1.second*fRNG.Gaus();
    fCurrentBeamShiftP2 = fBeamShiftCovP2.first + fBeamShiftCovP2.second*fRNG.Gaus();

    std::cout << gen_i << "Beam shift pars: " << fCurrentBeamShiftP0 <<
                 " " << fCurrentBeamShiftP1 << " " << fCurrentBeamShiftP2 <<
                 std::endl;
    */
    GenerateBeamShiftUniverse();

    for (auto it = new_samples.begin(); it != new_samples.end(); ++it) {
      for (size_t i = 0; i < it->second.size(); ++i) {
        for (size_t j = 0; j < it->second[i].size(); ++j) {
          it->second[i][j].Reset();
        }
      }
    }

    RefillMCSamples(events, new_samples, signal_sample_checks, beam_energy_bins,
                    signal_pars, flux_pars, syst_pars, g4rw_pars, fit_under_over,
                    tie_under_over, use_beam_inst_P);
    //To-do renormalize

    double total_varied_flux = 0.;
    double total_nominal_flux = 0.;
    std::vector<double> varied_flux(beam_energy_bins.size() - 1, 0.);
    std::vector<double> nominal_flux(beam_energy_bins.size() - 1, 0.);

    //Go through and determine the scaled flux
    //Two styles:
    //--Scaling to true incident momentum for pions/primary only
    //--Scaling to reco beam inst momentum for all types of particles
    //----(Set by use_beam_inst_P)

    for (auto it = fluxes_by_sample.begin(); it != fluxes_by_sample.end(); ++it) {
      int sample_ID = it->first;
      std::vector<std::vector<ThinSliceSample>> & samples_2D
          = new_samples[sample_ID];
      for (size_t i = 0; i < samples_2D.size(); ++i) {
        std::vector<ThinSliceSample> & samples = samples_2D[i];
        for (size_t j = 0; j < samples.size(); ++j) {
          ThinSliceSample & sample = samples.at(j);
          int flux_type = sample.GetFluxType();
          if (use_beam_inst_P || flux_type == 1) {
            nominal_flux[i] += fluxes_by_sample[sample_ID][i][j];
            varied_flux[i] += sample.GetVariedFlux();
          }
          total_nominal_flux += fluxes_by_sample[sample_ID][i][j];
          total_varied_flux += sample.GetVariedFlux();
        }
      }
    }

    double total_flux_factor = (total_nominal_flux/total_varied_flux);

    if (!use_beam_inst_P) {
      double nominal_primary = 0., varied_primary = 0.;
      for (size_t i = 0; i < nominal_flux.size(); ++i) {
        nominal_primary += nominal_flux[i];
        varied_primary += varied_flux[i];
      }

      std::vector<double> flux_factor = nominal_flux;
      for (size_t i = 0; i < flux_factor.size(); ++i) {
        flux_factor[i] = (nominal_flux[i] > 1.e-7 ?
                          flux_factor[i] / varied_flux[i] :
                          0.)*(varied_primary/nominal_primary);
      }

      for (auto it = new_samples.begin(); it != new_samples.end(); ++it) {
        for (size_t i = 0; i < it->second.size(); ++i) {
          for (size_t j = 0; j < it->second[i].size(); ++j) {
            int flux_type = it->second[i][j].GetFluxType();

            it->second[i][j].SetFactorAndScale(
                total_flux_factor*(flux_type == 1 ? flux_factor[i] : 1.));
          }
        }
      }
    }
    else {
      std::vector<double> flux_factor = nominal_flux;
      std::vector<double> new_var_flux(flux_factor.size());
      for (size_t i = 0; i < flux_factor.size(); ++i) {
        flux_factor[i] = (nominal_flux[i] > 1.e-7 ?
                          flux_factor[i] / varied_flux[i] :
                          0.);
      }

      for (auto it = new_samples.begin(); it != new_samples.end(); ++it) {
        for (size_t i = 0; i < it->second.size(); ++i) {
          for (size_t j = 0; j < it->second[i].size(); ++j) {
            it->second[i][j].SetFactorAndScale(flux_factor[i]);
            new_var_flux[i] += it->second[i][j].GetVariedFlux();
          }
        }
      }

    }

    //loop over varied/new samples and calculate bin_i_j as above
    //
    //
    //New samples -- get total histograms
    std::vector<double> new_vals;
    bool push_back = true;
    //double total_65 = 0.;
    for (auto it = new_samples.begin(); it != new_samples.end(); ++it) {
      for (size_t i = 0; i < it->second.size(); ++i) {
        for (size_t j = 0; j < it->second[i].size(); ++j) {
          const std::map<int, TH1*> & sel_hists
              = it->second[i][j].GetSelectionHists();
          int bin = 0;
          for (auto it2 = sel_hists.begin(); it2 != sel_hists.end(); ++it2) {
            if (it2->first == 6 || it2->first == 5) continue; 

            for (int k = 1; k <= it2->second->GetNbinsX(); ++k) {
              if (push_back) {
                new_vals.push_back(it2->second->GetBinContent(k));
              }
              else {
                new_vals[bin] += it2->second->GetBinContent(k);
              }
              ++bin;
            }
          }
          if (push_back) push_back = false;
        }
      }
    }


    all_new_vals.push_back(new_vals);
    for (size_t i = 0; i < new_vals.size(); ++i) {
      mean_new_vals[i] += new_vals[i];
    }
  }

  for (size_t i = 0; i < mean_new_vals.size(); ++i) {
    mean_new_vals[i] /= fNCovarianceGens;
  }

  std::vector<double> std_dev(mean_new_vals.size());
  for (auto & vals: all_new_vals) {
    for (size_t i = 0; i < vals.size(); ++i) {
      std_dev[i] += (vals[i] - mean_new_vals[i])*(vals[i] - mean_new_vals[i]);
    }
  }

  for (auto & sd : std_dev) {
    sd = sqrt(sd/fNCovarianceGens);
  }

  std::vector<std::vector<double>> good_vals;

  for (size_t gen_i = 0; gen_i < fNCovarianceGens; ++gen_i) {
    bool all_good = true;
    for (size_t i = 0; i < mean_new_vals.size(); ++i) {
      double nsigma = fabs((mean_new_vals[i] - all_new_vals[gen_i][i])/std_dev[i]);
      if (nsigma > 5.) {
        std::cout << gen_i << " " << i << " " << nsigma << std::endl;
        all_good = false;
        break;
      }
    }
    if (all_good)
      good_vals.push_back(all_new_vals[gen_i]);
  }

  for (size_t gen_i = 0; gen_i < good_vals.size(); ++gen_i) {
  //for (size_t gen_i = 0; gen_i < fNCovarianceGens; ++gen_i) {
    for (size_t i = 0; i < mean_new_vals.size(); ++i) {
      double val_i = (mean_new_vals[i] - good_vals[gen_i][i])/mean_new_vals[i];
      for (size_t j = 0; j < mean_new_vals.size(); ++j) {
        double val_j = (mean_new_vals[j] - good_vals[gen_i][j])/mean_new_vals[j];
        //std::cout << i << " " << j << " " << val_i << " " << val_j << std::endl;
        output_cov_mat[i][j] += val_i*val_j/fNCovarianceGens;
        output_cov.SetBinContent(i+1, j+1, (output_cov.GetBinContent(i+1, j+1) +
                                            val_i*val_j/fNCovarianceGens));
        //std::cout << "\t" << output_cov_mat[i][j] << " " <<
        //             output_cov.GetBinContent(i+1, j+1) << std::endl;
      }
    }
  }

  //Test decomposing it
  TDecompChol test_chol(output_cov_mat);
  bool decomposed = test_chol.Decompose();
  if (!decomposed) {
    std::cout << "WARNING: THE CONSTRUCTED COVARIANCE WAS NOT ABLE TO BE DECOMPOSED" << std::endl;
  }

  TH2D * output_corr = (TH2D*)output_cov.Clone("corr_hist");
  for (int i = 1; i <= output_corr->GetNbinsX(); ++i) {
    double val_i = output_cov.GetBinContent(i, i);
    for (int j = 1; j <= output_corr->GetNbinsX(); ++j) {
      double val_j = output_cov.GetBinContent(j, j);
      double content = output_corr->GetBinContent(i, j);
      output_corr->SetBinContent(
          i, j, content/sqrt(val_i*val_j));
    }
  }

  output_cov.Write();
  output_corr->Write();
  output_cov_mat.Write("cov_mat");

  TH1D centrals_hist("centrals_hist", "", mean_new_vals.size(), 0,
                     mean_new_vals.size());
  TVectorD centrals(mean_new_vals.size());
  for (size_t i = 0; i < mean_new_vals.size(); ++i) {
    centrals[i] = mean_new_vals[i]/nominal_vals[i];
    centrals_hist.SetBinContent(i+1, centrals[i]);
    centrals_hist.SetBinError(i+1, sqrt(output_cov_mat[i][i]));
  }
  centrals.Write("centrals");
  centrals_hist.Write();


  auto * dir = (TDirectory*)output_cov_file.mkdir("dists");
  dir->cd();
  TVectorD nom_out(nominal_vals.size());
  for (size_t i = 0; i < nominal_vals.size(); ++i) {
    nom_out(i) = nominal_vals[i]; 
  }
  nom_out.Write("dist_nom");

  int count = 0;
  for (auto & vals : all_new_vals) {
    TVectorD vals_out(vals.size());
    for (size_t i = 0; i < vals.size(); ++i) {
      vals_out(i) = vals[i]; 
    }

    std::string name = "dist_" + std::to_string(count);
    vals_out.Write(name.c_str());
    ++count;
  }

  count = 0;
  for (auto & vals : good_vals) {
    TVectorD vals_out(vals.size());
    for (size_t i = 0; i < vals.size(); ++i) {
      vals_out(i) = vals[i]; 
    }

    std::string name = "good_dist_" + std::to_string(count);
    vals_out.Write(name.c_str());
    ++count;
  }

  output_cov_file.Close();
  //Save the covariances in a new file
  //
  fBeamShiftCovRoutineActive = false;
}

double protoana::AbsCexDriver::GetBeamShiftDelta(const std::vector<double> & energies) {
  if (energies.size() == 0) return 0.;

  //Currently, Sungbin's work is in P_reco_BI - P_range_FF vs P_reco_BI
  //so we need to transform to P_reco_BI vs P_range_FF -- P_true_FF is approximately equal
  //
  //x - y = p2*x^2 + p1*x + p0
  //--> y = -p2*x^2 + (1 - p1)*x - p0
  //
  //Invert --> x = -p2*y^2 + (1 - p1)*y - p0
  //Find root of: 0 = -p2*y^2 + (1 - p1)*y - (p0 + x)
  //
  //Solve quad. eqn
  // --> y = ((1 - p1) +/- sqrt((1 - p1)^2 - 4*p2*(p0 + x)))/(2*p2)
  //
  // if p2 < 0., take positive soln., else take negative
  //
  // Check if the inverse is defined in the domain of P_true_FF
  // Don't generate a shift if not (set to 0.)
  // if p2 < 0, p0 > 0 --> defined for x > 0. --> all good
  //    p2 < 0, p0 < 0 --> defined for x >= (1 - p1)^2/(4*p2) - p0
  //    p2 > 0         --> defined for x <= (1 - p1)^2/(4*p2) - p0
  //
  // Also we have to subtract the nominal P_reco_BI vs P_true_FF

  double energy = energies[0] + 139.57;
  double true_P = sqrt(energy*energy - 139.57*139.57);

  if (fCurrentBeamShiftP2 < 0. && fCurrentBeamShiftP0 < 0.) {
    double domain_edge = ((1 - fCurrentBeamShiftP1)*(1 - fCurrentBeamShiftP1) /
                          (4*fCurrentBeamShiftP2)) - fCurrentBeamShiftP0;
    if (domain_edge > true_P) {
      std::cout << "P2 < 0" << std::endl;
      std::cout << "Skipping " << true_P << " " << domain_edge << std::endl;
      return 0.;
    }
  }
  else if (fCurrentBeamShiftP2 > 0.) {
    double domain_edge = ((1 - fCurrentBeamShiftP1)*(1 - fCurrentBeamShiftP1) /
                          (4*fCurrentBeamShiftP2)) - fCurrentBeamShiftP0;
    if (domain_edge < true_P) {
      std::cout << "P2 > 0" << std::endl;
      std::cout << "Skipping " << true_P << " " << domain_edge << std::endl;
      return 0.;
    }
  }

  double nominal_reco = fBeamShiftNominalP0 +
                        fBeamShiftNominalP1*true_P +
                        fBeamShiftNominalP2*true_P*true_P;

  double factor = (fCurrentBeamShiftP2 < 0. ? 1. : -1.); //If p2 < 0. use positive soln.
  double varied_reco
      = ((1. - fCurrentBeamShiftP1) +
         factor*sqrt((1. - fCurrentBeamShiftP1)*
                     (1. - fCurrentBeamShiftP1) -
                     4.*fCurrentBeamShiftP2*(fCurrentBeamShiftP0 + true_P)))*
        (1./(2*fCurrentBeamShiftP2));
  //std::cout << true_P << " " << nominal_reco << " " << varied_reco << " " <<
  //             varied_reco - nominal_reco << std::endl;
  //double shift = fCurrentBeamShiftP0 +
  //               fCurrentBeamShiftP1*energies[0] +
  //               fCurrentBeamShiftP2*energies[0]*energies[0];
  return 1.e-3*(varied_reco - nominal_reco);
  //return shift*1.e-3;
}

void protoana::AbsCexDriver::GenerateBeamShiftUniverse() {

    if (fBeamShiftCovUseInput) {
      TVectorD rand(3);
      for (size_t i = 0; i < 3; ++i) rand[i] = fRNG.Gaus();

      std::cout << "Using input covariance" << std::endl;
      for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
          //std::cout << fBeamShiftInputChol.GetU()[i][j] << " ";
          std::cout << (*fBeamShiftInputCovL)[i][j] << " ";
        }
        std::cout << std::endl;
      }
      //TVectorD rand_times_chol = fBeamShiftInputChol.GetU()*rand;
      TVectorD rand_times_chol = (*fBeamShiftInputCovL)*rand;
      for (size_t i = 0; i < 3; ++i) {
        std::cout << "\t" << rand_times_chol[i] << std::endl;
      }
      fCurrentBeamShiftP0 =  (*fBeamShiftInputCentrals)[0] + rand_times_chol[0];
      fCurrentBeamShiftP1 =  (*fBeamShiftInputCentrals)[1] + rand_times_chol[1];
      fCurrentBeamShiftP2 =  (*fBeamShiftInputCentrals)[2] + rand_times_chol[2];
    }

    else {
      fCurrentBeamShiftP0 = fBeamShiftCovP0.first + fBeamShiftCovP0.second*fRNG.Gaus();
      fCurrentBeamShiftP1 = fBeamShiftCovP1.first + fBeamShiftCovP1.second*fRNG.Gaus();
      fCurrentBeamShiftP2 = fBeamShiftCovP2.first + fBeamShiftCovP2.second*fRNG.Gaus();
    }

    std::cout << /*gen_i <<*/ "Beam shift pars: " << fCurrentBeamShiftP0 <<
                 " " << fCurrentBeamShiftP1 << " " << fCurrentBeamShiftP2 <<
                 std::endl;
}

void protoana::AbsCexDriver::FillDataHistsFromSamples(
    ThinSliceDataSet & data_set,
    const std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    double & flux, std::vector<double> & fluxes_by_beam, bool fluctuate) {

  auto & selection_hists = data_set.GetSelectionHists();

  flux = 0.;
  //std::vector<double> temp_flux(/*fluxes_by_beam.size()*/);
  fluxes_by_beam.clear();
  for (auto it = selection_hists.begin(); it != selection_hists.end(); ++it) {
    it->second->Reset();
  }


  if (fluctuate) {
    //Get total
    double total = 0.;
    int a = 0;
    for (auto & sample : samples ) {
      for (size_t i = 0; i < sample.second.size(); ++i) {
        for (size_t j = 0; j < sample.second[i].size(); ++j) {
          const auto & hists = sample.second[i][j].GetSelectionHists();
          for (auto & hist : hists) {
            total += hist.second->Integral();
            std::cout << a << " " << i << " "  << j << " " << hist.second->Integral() << std::endl;
            //std::cout << temp_flux.size() << std::endl;
            //temp_flux[i] += hist.second->Integral();
          }
        }
      }
      ++a;
    }
    std::cout << "Total: " << total << std::endl;
    //for (auto & f : temp_flux) {std::cout << f << " ";} std::cout << std::endl;


    //First, iterate through all samples/hists. Vary the contents
    //Skip accordingly. Get varied total
    double varied_total = 0.;
    double nominal_skip_total = 0.;
    //std::map<int, double> nominal_skip_vals;
    //for (auto i : fToSkip) nominal_skip_vals[i] = 0.;

    for (auto & sample : samples ) {
      auto & sample_vec = sample.second;
      for (size_t i = 0; i < sample_vec.size(); ++i) {
        for (size_t j = 0; j < sample_vec[i].size(); ++j) {
          const auto & hists = sample_vec[i][j].GetSelectionHists();
          for (auto & hist : hists) {
            if (std::find(fToSkip.begin(), fToSkip.end(), hist.first)
                != fToSkip.end()) {
              nominal_skip_total += hist.second->Integral();
              //nominal_skip_vals[hist.first] += hist.second->Integral();
              continue;
            }
            for (int k = 1; k <= hist.second->GetNbinsX(); ++k) {
              double val = hist.second->GetBinContent(k);
              double rand = fRNG.PoissonD(val);
              hist.second->SetBinContent(k, rand);
            }
            varied_total += hist.second->Integral();
          }
        }
      }
    }

    double leftover = total - varied_total;
    double ratio = leftover/nominal_skip_total;
    std::cout << "Leftover: " << leftover << std::endl;
    std::cout << "Nominal skip total " << nominal_skip_total << std::endl;

    for (auto & sample : samples ) {
      auto & sample_vec = sample.second;
      for (size_t i = 0; i < sample_vec.size(); ++i) {
        for (size_t j = 0; j < sample_vec[i].size(); ++j) {
          const auto & hists = sample_vec[i][j].GetSelectionHists();
          for (auto k : fToSkip) {
            if (hists.find(k) == hists.end()) continue;
            //std::cout << "Scaling " << k << " " << hists.at(k)->Integral() <<
            //             " " << ratio << std::endl;
            hists.at(k)->Scale(ratio);
          }
        }
      }
    }
  }

  int a = 0;
  for (auto it = samples.begin(); it != samples.end(); ++it) {
    for (size_t i = 0; i < it->second.size(); ++i) {
      if (a == 0) fluxes_by_beam.push_back(0.);
      for (size_t j = 0; j < it->second[i].size(); ++j) {
        const auto & hists = it->second[i][j].GetSelectionHists();

        for (auto it2 = hists.begin(); it2 != hists.end(); ++it2) {

          selection_hists[it2->first]->Add(it2->second);
          flux += it2->second->Integral();
          fluxes_by_beam[i] += it2->second->Integral();
          //std::cout << "Adding " << it2->second->Integral() << " to " << i << std::endl;
        }
      }
    }
    ++a;
  }
  std::cout << "Flux: " << flux << std::endl;
}
