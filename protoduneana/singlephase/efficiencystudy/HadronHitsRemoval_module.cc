////////////////////////////////////////////////////////////////////////
// Class:       HadronHitsRemoval
// Plugin Type: producer (art v3_06_03)
// File:        HadronHitsRemoval_module.cc
//
// Generated at Tue Jul 13 22:36:12 2021 by Tingjun Yang using cetskelgen
// from cetlib version v3_11_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

#include "art_root_io/TFileService.h"
#include "TTree.h"
#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "TMath.h"
#include "TVector3.h"
using namespace std;
namespace pdsp {
  class HadronHitsRemoval;
}


class pdsp::HadronHitsRemoval : public art::EDProducer {
public:
  explicit HadronHitsRemoval(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  HadronHitsRemoval(HadronHitsRemoval const&) = delete;
  HadronHitsRemoval(HadronHitsRemoval&&) = delete;
  HadronHitsRemoval& operator=(HadronHitsRemoval const&) = delete;
  HadronHitsRemoval& operator=(HadronHitsRemoval&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;
  void beginJob() override;
  void endJob() override;
  void reset();
  
  bool PassBeamQualityCut() const;

private:
  // Declare member data here.
  TTree *fTree;
  double reco_beam_len_sce;
  int selected_track;
  protoana::ProtoDUNETrackUtils trackUtil;
  protoana::ProtoDUNEPFParticleUtils pfpUtil;
  // for beam quality cut of data
  const double beam_startX_data = -27.911;
  const double beam_startY_data = 424.364;
  const double beam_startZ_data = 3.77836;
  const double beam_startX_rms_data = 4.71128;
  const double beam_startY_rms_data = 5.16472;
  const double beam_startZ_rms_data = 1.10265;
  const double beam_angleX_data = 100.454;
  const double beam_angleY_data = 103.523;
  const double beam_angleZ_data = 17.8288;

  double reco_beam_calo_startX, reco_beam_calo_startY, reco_beam_calo_startZ;
  double reco_beam_calo_endX, reco_beam_calo_endY, reco_beam_calo_endZ;
  double beam_dx, beam_dy, beam_dz, beam_dxy, beam_costh;
};


pdsp::HadronHitsRemoval::HadronHitsRemoval(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
  // More initializers here.
{
  cout<<"$$$HadronHitsRemoval"<<endl;
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void pdsp::HadronHitsRemoval::beginJob(){
  cout<<"$$$beginJob"<<endl;
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("beamana","beam analysis tree");
  fTree->Branch("reco_beam_len_sce", &reco_beam_len_sce);
  fTree->Branch("beam_dx", &beam_dx);
  fTree->Branch("beam_dy", &beam_dy);
  fTree->Branch("beam_dz", &beam_dz);
  fTree->Branch("beam_dxy", &beam_dxy);
  fTree->Branch("beam_costh", &beam_costh);
  fTree->Branch("selected_track", &selected_track);
}

void pdsp::HadronHitsRemoval::produce(art::Event& evt)
{
  reset();
  string fPFParticleTag = "pandora";
  string fCalorimetryTagSCE = "pandoracalo";
  string fTrackerTag = "pandoraTrack";
  // Implementation of required member function here.
  // Add code to select beam tracks using Pandora information
  
  const std::map<unsigned int, std::vector<const recob::PFParticle*>> sliceMap
      = pfpUtil.GetPFParticleSliceMap(evt, fPFParticleTag);
  std::vector<std::vector<const recob::PFParticle*>> beam_slices;
  for(auto slice : sliceMap){
    for(auto particle : slice.second){
      bool added = false;
      if(pfpUtil.IsBeamParticle(*particle,evt, fPFParticleTag)){
        if (!added) {
          beam_slices.push_back(slice.second);
          added = true;
        }
      }
      if (!added) {continue;}
      //else {
      //  std::cout << "Not beam particle" << std::endl;
      //}
    }
  }
  if (beam_slices.size() > 0){
    const recob::PFParticle* particle = beam_slices.at(0).at(0);
    
    const recob::Track* thisTrack = pfpUtil.GetPFParticleTrack(*particle,evt,fPFParticleTag,fTrackerTag);

    //Primary Track Calorimetry
    if( thisTrack ){
      auto calo = trackUtil.GetRecoTrackCalorimetry(*thisTrack, evt, fTrackerTag,
                                                    fCalorimetryTagSCE);
      bool found_calo = false;
      size_t index = 0;
      for ( index = 0; index < calo.size(); ++index) {
        if (calo[index].PlaneID().Plane == 2) {
          found_calo = true;
          break;
        }
      }
      if (found_calo) {
        reco_beam_len_sce = calo[index].Range();
        auto theXYZPoints = calo[index].XYZ();
        std::sort(theXYZPoints.begin(), theXYZPoints.end(), [](auto a, auto b)
            {return (a.Z() < b.Z());}); // sort according to Z position
        if (theXYZPoints.size()) {
          reco_beam_calo_startX = theXYZPoints[0].X();
          reco_beam_calo_startY = theXYZPoints[0].Y();
          reco_beam_calo_startZ = theXYZPoints[0].Z();
          reco_beam_calo_endX = theXYZPoints.back().X();
          reco_beam_calo_endY = theXYZPoints.back().Y();
          reco_beam_calo_endZ = theXYZPoints.back().Z();
        }
        
        // beam quality requirement (for data)
        beam_dx = (reco_beam_calo_startX - beam_startX_data)/beam_startX_rms_data;
        beam_dy = (reco_beam_calo_startY - beam_startY_data)/beam_startY_rms_data;
        beam_dz = (reco_beam_calo_startZ - beam_startZ_data)/beam_startZ_rms_data;
        beam_dxy = sqrt(pow(beam_dx,2) + pow(beam_dy,2));
        
        TVector3 beamdir_data(cos(beam_angleX_data*TMath::Pi()/180),
                              cos(beam_angleY_data*TMath::Pi()/180),
                              cos(beam_angleZ_data*TMath::Pi()/180));
        beamdir_data = beamdir_data.Unit();
        
        TVector3 pt0(reco_beam_calo_startX,
                     reco_beam_calo_startY,
                     reco_beam_calo_startZ);
        TVector3 pt1(reco_beam_calo_endX,
                     reco_beam_calo_endY,
                     reco_beam_calo_endZ);
        TVector3 dir = pt1 - pt0;
        dir = dir.Unit();
        beam_costh = dir.Dot(beamdir_data);

        if (reco_beam_len_sce>100 && PassBeamQualityCut()){
          selected_track = 1;
          if (reco_beam_len_sce>230)
            selected_track = 2;
        }
      }
    }
  }

  cout<<"$$$reco_beam_len_sce "<<reco_beam_len_sce<<endl;
  fTree->Fill();
}

void pdsp::HadronHitsRemoval::endJob(){
  cout<<"$$$endJob"<<endl;
}

void pdsp::HadronHitsRemoval::reset(){
  cout<<"$$$reset"<<endl;
  reco_beam_len_sce = -999;
  beam_dx = -999;
  beam_dy = -999;
  beam_dz = -999;
  beam_dxy = -999;
  beam_costh = -999;
  reco_beam_calo_startX = -999;
  reco_beam_calo_startY = -999;
  reco_beam_calo_startZ = -999;
  reco_beam_calo_endX = -999;
  reco_beam_calo_endY = -999;
  reco_beam_calo_endZ = -999;
  selected_track = 0;
}

bool pdsp::HadronHitsRemoval::PassBeamQualityCut() const{ // cut on beam entrance location and beam angle

  if (beam_dxy<-1) return false;
  if (beam_dxy>3) return false;
  
  if (beam_dz<-3) return false;
  if (beam_dz>3) return false;
  
  if (beam_costh<0.95) return false;
  if (beam_costh>2) return false;
  
  return true;
}

DEFINE_ART_MODULE(pdsp::HadronHitsRemoval)
