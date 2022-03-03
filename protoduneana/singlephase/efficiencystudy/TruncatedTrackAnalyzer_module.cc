////////////////////////////////////////////////////////////////////////
// Class:       TruncatedTrackAnalyzer
// Plugin Type: analyzer (Unknown Unknown)
// File:        TruncatedTrackAnalyzer_module.cc
//
// Generated at Thu Mar  3 11:53:02 2022 by Jacob Calcutt using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "art_root_io/TFileService.h"
#include "TTree.h"

namespace pdsp {
  class TruncatedTrackAnalyzer;
}


class pdsp::TruncatedTrackAnalyzer : public art::EDAnalyzer {
public:
  explicit TruncatedTrackAnalyzer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TruncatedTrackAnalyzer(TruncatedTrackAnalyzer const&) = delete;
  TruncatedTrackAnalyzer(TruncatedTrackAnalyzer&&) = delete;
  TruncatedTrackAnalyzer& operator=(TruncatedTrackAnalyzer const&) = delete;
  TruncatedTrackAnalyzer& operator=(TruncatedTrackAnalyzer&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;
  void beginJob() override;

  void reset();
  void CheckEff(const art::Event & evt, int rightTrackID);

private:
  protoana::ProtoDUNEPFParticleUtils fPFPUtil;

  //Output
  TTree * fTree;
  int reco_beam_type, reco_beam_trackID;
  int for_truncation_method;
  int run, subrun, event;
  int has_truncated_track;
  int n_truncated_hits;
  std::vector<double> truncated_hits_X, truncated_hits_Y, truncated_hits_Z;

  //Configuration
  std::string fPFParticleTag, fTrackerTag, fShowerTag, fHitModuleTag;

};


pdsp::TruncatedTrackAnalyzer::TruncatedTrackAnalyzer(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fPFParticleTag(p.get<std::string>("PFParticleTag")),
    fTrackerTag(p.get<std::string>("TrackerTag")),
    fShowerTag(p.get<std::string>("ShowerTag")) {
}

void pdsp::TruncatedTrackAnalyzer::analyze(art::Event const& evt) {
  reset();

  run = evt.run();
  subrun = evt.subRun();
  event = evt.id().event();

  const std::vector<const recob::PFParticle*> beam_pfps
      = fPFPUtil.GetPFParticlesFromBeamSlice(evt, fPFParticleTag);
  std::cout << "Got " << beam_pfps.size() << " beam pfps" << std::endl;
  if (beam_pfps.size()) {
    const recob::Track* thisTrack = fPFPUtil.GetPFParticleTrack(
        *(beam_pfps[0]), evt, fPFParticleTag, fTrackerTag);
    const recob::Shower * thisShower = fPFPUtil.GetPFParticleShower(
        *(beam_pfps[0]), evt, fPFParticleTag, fShowerTag);
    
    if (thisTrack) {
      reco_beam_type = 13;
      reco_beam_trackID = thisTrack->ID();
    }
    else if (thisShower) {
      reco_beam_type = 11;
    }
    else {
      reco_beam_type = -999;
    }

  }
  CheckEff(evt, reco_beam_trackID);

  fTree->Fill();
}

void pdsp::TruncatedTrackAnalyzer::reset() {
  reco_beam_type = -999;
  event = -999;
  run = -999;
  subrun = -999;
  reco_beam_trackID = -999;
  for_truncation_method = -999;
  has_truncated_track = 0;
  n_truncated_hits = 0;
  truncated_hits_X.clear();
  truncated_hits_Y.clear();
  truncated_hits_Z.clear();
}

void pdsp::TruncatedTrackAnalyzer::beginJob() {
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("tree", "");

  fTree->Branch("run", &run);
  fTree->Branch("subrun", &subrun);
  fTree->Branch("event", &event);
  fTree->Branch("reco_beam_type", &reco_beam_type);
  fTree->Branch("reco_beam_trackID", &reco_beam_trackID);
  fTree->Branch("for_truncation_method", &for_truncation_method);
  fTree->Branch("has_truncated_track", &has_truncated_track);
  fTree->Branch("n_truncated_hits", &n_truncated_hits);
  fTree->Branch("truncated_hits_X", &truncated_hits_X);
  fTree->Branch("truncated_hits_Y", &truncated_hits_Y);
  fTree->Branch("truncated_hits_Z", &truncated_hits_Z);
}

void pdsp::TruncatedTrackAnalyzer::CheckEff(const art::Event & evt, int rightTrackID) {
  auto hitHandler = evt.getValidHandle<std::vector<recob::Hit> >("hitpdune");
  auto spHandler = evt.getValidHandle< std::vector<recob::SpacePoint> >("hitpdune");
  art::FindManyP<recob::Hit> hitFromSP(spHandler, evt, "hitpdune");
  art::FindManyP<recob::Track> trackFromHit(hitHandler, evt, "pandoraTrack");

  art::FindManyP<recob::SpacePoint> pandoraSPFromHit(hitHandler, evt, "pandora");

  std::cout << "hitsFromSP: " << hitFromSP.size() << std::endl;
  if (hitFromSP.size() > 0) {
    has_truncated_track = 1;
    for_truncation_method = 0;
    auto const & hits = hitFromSP.at(0);
    n_truncated_hits = hits.size();
    std::vector<int> keys;

    //Get the calorimetry object out


    for (auto const & hit: hits){
      auto const & tracks = trackFromHit.at(hit.key());
      if (!tracks.empty()){
        keys.push_back(tracks[0].key());
      }

      auto const & sps = pandoraSPFromHit.at(hit.key());
      for (auto sp : sps) {
        truncated_hits_X.push_back(sp->XYZ()[0]);
        truncated_hits_Y.push_back(sp->XYZ()[1]);
        truncated_hits_Z.push_back(sp->XYZ()[2]);
      }
    }
    
    // find the most common element (mode) in keys
    int repetition = 0;
    int mode = -2;
    std::map<int,int> mmap;
    for (auto vi: keys) {
      mmap[vi]++;
      if (mmap[vi] > repetition) {
        repetition = mmap[vi];
        mode = vi;
      }
    }
    
    if (repetition > hits.size()/2.) { // if repetition of mode in keys > 1/2 size of tagged hits
      for_truncation_method = 1; // reconstructed
      if (mode == rightTrackID) {
        for_truncation_method = 2; // identified
      }
    }
  }
}

DEFINE_ART_MODULE(pdsp::TruncatedTrackAnalyzer)
