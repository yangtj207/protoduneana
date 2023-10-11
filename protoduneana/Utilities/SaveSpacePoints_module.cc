////////////////////////////////////////////////////////////////////////
// Class:       SaveSpacePoints
// Plugin Type: analyzer (art v2_11_02)
// File:        SaveSpacePoints_module.cc
//
// Generated at Wed Jul 11 09:18:07 2018 by Tingjun Yang using cetskelgen
// from cetlib version v3_03_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PointCharge.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RawData/RDTimeStamp.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "dunecore/DuneObj/ProtoDUNEBeamEvent.h"
#include "protoduneana/Utilities/ProtoDUNEBeamlineUtils.h"

// ROOT includes
#include "TTree.h"
#include "TTimeStamp.h"

namespace proto {
  class SaveSpacePoints;
}


class proto::SaveSpacePoints : public art::EDAnalyzer {
public:
  explicit SaveSpacePoints(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SaveSpacePoints(SaveSpacePoints const &) = delete;
  SaveSpacePoints(SaveSpacePoints &&) = delete;
  SaveSpacePoints & operator = (SaveSpacePoints const &) = delete;
  SaveSpacePoints & operator = (SaveSpacePoints &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;

private:

  const art::InputTag fSpacePointModuleLabel;
  const art::InputTag fBeamModuleLabel;
  const art::InputTag fTrackModuleLabel;
  const art::InputTag fTimeDecoderModuleLabel;

  protoana::ProtoDUNEBeamlineUtils fBeamlineUtils;

  TTree *fTree;
  // Run information
  int run;
  int subrun;
  int event;
  int trigger;
  double evttime;

  // space point information
  std::vector<double> vx;
  std::vector<double> vy;
  std::vector<double> vz;
  std::vector<double> vcharge;
  std::vector<int> vtpc;
  std::vector<int> vtrackid;
  std::vector<int> vpdg;
  std::vector<int> vg4id;
  std::vector<int> vorigin;

  // beam information
  std::vector<double> beamPosx;
  std::vector<double> beamPosy;
  std::vector<double> beamPosz;
  
  std::vector<double> beamDirx;
  std::vector<double> beamDiry;
  std::vector<double> beamDirz;

  std::vector<double> beamMomentum;

  double tof;
  short ckov0status;
  short ckov1status;

};


proto::SaveSpacePoints::SaveSpacePoints(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  fSpacePointModuleLabel(p.get< art::InputTag >("SpacePointModuleLabel")),
  fBeamModuleLabel(p.get< art::InputTag >("BeamModuleLabel")),
  fTrackModuleLabel(p.get< art::InputTag >("TrackModuleLabel")),
  fTimeDecoderModuleLabel(p.get< art::InputTag >("TimeDecoderModuleLabel")),
  fBeamlineUtils(p.get<fhicl::ParameterSet>("BeamlineUtils"))
{}

void proto::SaveSpacePoints::analyze(art::Event const & evt)
{
  // Implementation of required member function here.
  run = evt.run();
  subrun = evt.subRun();
  event = evt.id().event();
  art::Timestamp ts = evt.time();
  //std::cout<<ts.timeHigh()<<" "<<ts.timeLow()<<std::endl;
  if (ts.timeHigh() == 0){
    TTimeStamp tts(ts.timeLow());
    evttime = tts.AsDouble();
  }
  else{
    TTimeStamp tts(ts.timeHigh(), ts.timeLow());
    evttime = tts.AsDouble();
  }
  vx.clear();
  vy.clear();
  vz.clear();
  vcharge.clear();
  vtpc.clear();
  vtrackid.clear();
  vpdg.clear();
  vg4id.clear();
  vorigin.clear();
  beamPosx.clear();
  beamPosy.clear();
  beamPosz.clear();
  beamDirx.clear();
  beamDiry.clear();
  beamDirz.clear();
  beamMomentum.clear();

  if (evt.isRealData()){
    if (fBeamlineUtils.IsGoodBeamlineTrigger(evt)){

      //Access the Beam Event
      auto beamHandle = evt.getValidHandle<std::vector<beam::ProtoDUNEBeamEvent>>("beamevent");
      
      std::vector<art::Ptr<beam::ProtoDUNEBeamEvent>> beamVec;
      if( beamHandle.isValid()){
        art::fill_ptr_vector(beamVec, beamHandle);
      }
      
      const beam::ProtoDUNEBeamEvent & beamEvent = *(beamVec.at(0)); //Should just have one
      
      //Access momentum
      const std::vector< double > & the_momenta = beamEvent.GetRecoBeamMomenta();
      std::cout << "Number of reconstructed momenta: " << the_momenta.size() << std::endl;
      
      beamMomentum.insert( beamMomentum.end(), the_momenta.begin(), the_momenta.end() );
      
      tof = -1;
      ckov0status = -1;
      ckov1status = -1;
      
      //Access time of flight
      const std::vector< double > & the_tofs  = beamEvent.GetTOFs();
      
      if( the_tofs.size() > 0){
        tof = the_tofs[0];
      }
      ckov0status = beamEvent.GetCKov0Status();
      ckov1status = beamEvent.GetCKov1Status();
    }
  
    trigger = -1;
    art::ValidHandle<std::vector<raw::RDTimeStamp>> timeStamps = evt.getValidHandle<std::vector<raw::RDTimeStamp>>(fTimeDecoderModuleLabel);

    // Check that we have good information
    if(timeStamps.isValid() && timeStamps->size() == 1){
      // Access the trigger information. Beam trigger flag = 0xc
      const raw::RDTimeStamp& timeStamp = timeStamps->at(0);
      trigger = timeStamp.GetFlags();
    }
  }

  std::vector< art::Ptr<recob::SpacePoint> > sps;
  auto spsHandle = evt.getHandle< std::vector<recob::SpacePoint> >(fSpacePointModuleLabel);
  if (spsHandle) art::fill_ptr_vector(sps, spsHandle);

  std::vector< art::Ptr<recob::PointCharge> > pcs;
  auto pcsHandle = evt.getHandle< std::vector<recob::PointCharge> >(fSpacePointModuleLabel);
  if (pcsHandle) art::fill_ptr_vector(pcs, pcsHandle);

  art::FindManyP<recob::Hit> fmhsp(spsHandle, evt, fSpacePointModuleLabel);

  //Services
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);

  for (size_t i = 0; i<sps.size(); ++i){
    vx.push_back(sps[i]->XYZ()[0]);
    vy.push_back(sps[i]->XYZ()[1]);
    vz.push_back(sps[i]->XYZ()[2]);
    vcharge.push_back(pcs[i]->charge());
    vtrackid.push_back(-1);
    auto const& hits = fmhsp.at(i);
    if (!evt.isRealData()){
      int TrackID = 0;
      std::map<int,double> trkide;
      for (auto const & hit : hits){
        std::vector<sim::TrackIDE> TrackIDs = bt_serv->HitToTrackIDEs(clockData, hit);
        for(size_t e = 0; e < TrackIDs.size(); ++e){
          trkide[TrackIDs[e].trackID] += TrackIDs[e].energy;
        }
      }
      // Work out which IDE despoited the most charge in the hit if there was more than one.
      double maxe = -1;
      // double tote = 0; // unused
      for (std::map<int,double>::iterator ii = trkide.begin(); ii!=trkide.end(); ++ii){
        // tote += ii->second; // unused
        if ((ii->second)>maxe){
          maxe = ii->second;
          TrackID = ii->first;
        }
      }
      // Now have trackID, so get PdG code and T0 etc.
      vg4id.push_back(TrackID);
      const simb::MCParticle *particle = pi_serv->TrackIdToParticle_P(TrackID);
      if (particle){
        vpdg.push_back(particle->PdgCode());
        vorigin.push_back(pi_serv->ParticleToMCTruth_P(particle)->Origin());
      }
      else{
        vpdg.push_back(0);
        vorigin.push_back(0);
      }
    }
    else{
      vg4id.push_back(0);
      vpdg.push_back(0);
      vorigin.push_back(0);
    }
    int spTPC = -1;
    for (auto const & hit : hits){
      spTPC = hit->WireID().TPC;
    }
    vtpc.push_back(spTPC);
  }

  std::vector< art::Ptr<recob::Track> > trks;
  auto trkHandle = evt.getHandle< std::vector<recob::Track> >(fTrackModuleLabel);
  if (trkHandle) art::fill_ptr_vector(trks, trkHandle);

  for (size_t i = 0; i<trks.size(); ++i){
    auto & trk = trks[i];
    for (size_t j = 0; j<trk->NPoints(); ++j){
      if (trk->HasValidPoint(j)){
        vx.push_back(trk->TrajectoryPoint(j).position.X());
        vy.push_back(trk->TrajectoryPoint(j).position.Y());
        vz.push_back(trk->TrajectoryPoint(j).position.Z());
        vcharge.push_back(0);
        vtpc.push_back(-1);
        vtrackid.push_back(trk->ID());
        vpdg.push_back(-1);
        vg4id.push_back(-1);
        vorigin.push_back(-1);
      }
    }
  }

  fTree->Fill();
}

void proto::SaveSpacePoints::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("spt","space point tree");
  fTree->Branch("run",&run,"run/I");
  fTree->Branch("subrun",&subrun,"subrun/I");
  fTree->Branch("event",&event,"event/I");
  fTree->Branch("trigger",&trigger,"trigger/I");
  fTree->Branch("evttime",&evttime,"evttime/D");
  fTree->Branch("vx",&vx);
  fTree->Branch("vy",&vy);
  fTree->Branch("vz",&vz);
  fTree->Branch("vcharge",&vcharge);
  fTree->Branch("vtpc",&vtpc);
  fTree->Branch("vtrackid",&vtrackid);
  fTree->Branch("vpdg",&vpdg);
  fTree->Branch("vg4id",&vg4id);
  fTree->Branch("vorigin",&vorigin);
  fTree->Branch("beamPosx",&beamPosx);
  fTree->Branch("beamPosy",&beamPosy);
  fTree->Branch("beamPosz",&beamPosz);
  fTree->Branch("beamDirx",&beamDirx);
  fTree->Branch("beamDiry",&beamDiry);
  fTree->Branch("beamDirz",&beamDirz);
  fTree->Branch("beamMomentum",&beamMomentum);
  fTree->Branch("tof", &tof, "tof/D");
  fTree->Branch("ckov0status", &ckov0status, "ckov0status/S");
  fTree->Branch("ckov1status", &ckov1status, "ckov1status/S");
}

DEFINE_ART_MODULE(proto::SaveSpacePoints)
