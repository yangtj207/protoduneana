////////////////////////////////////////////////////////////////////////
// Class:       BoxBeamHitsRemoval
// Plugin Type: producer (Unknown Unknown)
// File:        BoxBeamHitsRemoval_module.cc
//
// Generated at Wed Mar  2 12:11:06 2022 by Jacob Calcutt using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardata/ArtDataHelper/HitCreator.h"

#include <memory>

namespace pdsp {
  class BoxBeamHitsRemoval;
}


class pdsp::BoxBeamHitsRemoval : public art::EDProducer {
public:
  explicit BoxBeamHitsRemoval(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  BoxBeamHitsRemoval(BoxBeamHitsRemoval const&) = delete;
  BoxBeamHitsRemoval(BoxBeamHitsRemoval&&) = delete;
  BoxBeamHitsRemoval& operator=(BoxBeamHitsRemoval const&) = delete;
  BoxBeamHitsRemoval& operator=(BoxBeamHitsRemoval&&) = delete;

  void produce(art::Event& e) override;

private:
  std::vector<std::pair<double, double>> fRegionsY;
  double fUpperZCut;
  bool fFlipTracks;
  std::string fPFParticleTag, fTrackerTag, fShowerTag, fHitModuleTag;
  protoana::ProtoDUNEPFParticleUtils fPFPUtil;
};


pdsp::BoxBeamHitsRemoval::BoxBeamHitsRemoval(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fRegionsY(p.get<std::vector<std::pair<double, double>>>("RegionsY")),
    fUpperZCut(p.get<double>("UpperZCut")),
    fFlipTracks(p.get<bool>("FlipTracks")),
    fPFParticleTag(p.get<std::string>("PFParticleTag")),
    fTrackerTag(p.get<std::string>("TrackerTag")),
    fShowerTag(p.get<std::string>("ShowerTag")),
    fHitModuleTag(p.get<std::string>("HitModuleTag", "hitpdune")) {
  recob::HitCollectionCreator::declare_products(producesCollector(), "", true, false);
  produces<std::vector<recob::SpacePoint>>();
  produces<art::Assns<recob::Hit, recob::SpacePoint>>(); // this space point will only be used as a tag
}

void pdsp::BoxBeamHitsRemoval::produce(art::Event& evt) {
  recob::HitCollectionCreator hcol(evt, true, false);
  std::unique_ptr<std::vector<recob::SpacePoint>> spvp(new std::vector<recob::SpacePoint>);
  auto assns = std::make_unique<art::Assns<recob::Hit, recob::SpacePoint>>();

  //Create vector to store hits that will eventually be removed
  std::vector<art::Ptr<recob::Hit>> hits_to_remove;


  //Loop over all PFPs in event. Look for those whose start falls in the box beam region
  auto pfpVec = evt.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleTag);
  for (const recob::PFParticle & pfp : (*pfpVec)) {
    const recob::Track* thisTrack = fPFPUtil.GetPFParticleTrack(
        pfp, evt, fPFParticleTag, fTrackerTag);
    const recob::Shower * thisShower = fPFPUtil.GetPFParticleShower(
        pfp, evt, fPFParticleTag, fShowerTag);

    if (!thisTrack && !thisShower) continue;

    bool near_box_beam = false;
 
    if (thisTrack) {
      double startY = thisTrack->Trajectory().Start().Y();
      double startZ = thisTrack->Trajectory().Start().Z();
      double endY = thisTrack->Trajectory().End().Y();
      double endZ = thisTrack->Trajectory().End().Z();

      //if flipping tracks, check if endZ > startZ
      if (fFlipTracks && (startZ > endZ)) {
        startZ = endZ;
        startY = endY;
      }
      
      if (startZ < fUpperZCut) {
        for (auto region : fRegionsY) {
          if (region.first < startY && startY < region.second) {
            near_box_beam = true;
            break;
          }
        }
      }
    }

    if (thisShower) {
      double startY = thisShower->ShowerStart().Y();  
      double startZ = thisShower->ShowerStart().Z();
      if (startZ < fUpperZCut) {
        for (auto region : fRegionsY) {
          if (region.first < startY && startY < region.second) {
            near_box_beam = true;
            break;
          }
        }
      }
    }

    if (!near_box_beam) continue;

    const std::vector<art::Ptr<recob::Hit>> pfp_hits = fPFPUtil.GetPFParticleHits_Ptrs(
        pfp, evt, fPFParticleTag);
    for (auto hit : pfp_hits) {
      hits_to_remove.push_back(hit);
    }
  }

  //Get all hits in the event
  auto hitsHandle = evt.getValidHandle<std::vector<recob::Hit>>(fHitModuleTag);
  art::FindOneP<recob::Wire> channelHitWires(hitsHandle, evt, fHitModuleTag);
  std::vector<art::Ptr<recob::Hit>> eventHits;
  art::fill_ptr_vector(eventHits, hitsHandle);

  // fill hits
  for (size_t i = 0; i < eventHits.size(); ++i){
    auto hit = eventHits[i];
    auto it = find(hits_to_remove.begin(), hits_to_remove.end(), hit);
    if (it == hits_to_remove.end()) { // not be to removed
      hcol.emplace_back(std::move(*hit), channelHitWires.at(i));
    }
  }

  hcol.put_into(evt);
  evt.put(std::move(spvp));
  evt.put(std::move(assns));
}

DEFINE_ART_MODULE(pdsp::BoxBeamHitsRemoval)
