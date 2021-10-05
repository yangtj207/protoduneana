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
#include <memory>

#include <string>
#include <vector>
#include <utility>
#include <memory>
#include <iostream>

// Framework includes
#include "canvas/Utilities/InputTag.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/EDProducer.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft Includes
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardata/ArtDataHelper/HitCreator.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

namespace pdsp {
  class HadronHitsRemoval;
}
using namespace std;


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
  geo::GeometryCore const* fGeom;
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
  // Declare member data here.
};


pdsp::HadronHitsRemoval::HadronHitsRemoval(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
  // More initializers here.
{
  cout<<"$$$HadronHitsRemoval$$$"<<endl;
  // Call appropriate produces<>() functions here.
  recob::HitCollectionCreator::declare_products(producesCollector(), "", true, false);
  fGeom = &*(art::ServiceHandle<geo::Geometry>());
  //produces<art::Assns<recob::Hit, recob::SpacePoint>>();
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
  cout<<"#####EvtNo."<<evt.id().event()<<endl;
  reset();
  string fPFParticleTag = "pandora";
  string fCalorimetryTagSCE = "pandoracalo";
  string fTrackerTag = "pandoraTrack";
  // Implementation of required member function here.
  // Add code to select beam tracks using Pandora information
  
  recob::HitCollectionCreator hcol(evt, true, false);
  
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
          if (reco_beam_len_sce>230 && reco_beam_len_sce<500)
            selected_track = 2;
        }
      }
      
      string fHitModuleLabel = "hitpdune";
      string fSpModuleLabel = "reco3d";
      
      double sel_len = 50.; // shorten to how long (cm)
      
      if (selected_track == 2) { // selected long tracks (not smaller than 5 m)
        auto hitsHandle = evt.getValidHandle< std::vector<recob::Hit> >(fHitModuleLabel);
        auto spHandle = evt.getValidHandle< std::vector<recob::SpacePoint> >(fSpModuleLabel);
        art::FindOneP<recob::Wire>   channelHitWires    (hitsHandle, evt, fHitModuleLabel);
        art::FindManyP< recob::SpacePoint > spFromHit(hitsHandle, evt, fSpModuleLabel);
        //recob::HitCollectionCreator hcol(evt, true, true);
        std::vector< art::Ptr<recob::Hit> > eventHits;
        art::fill_ptr_vector(eventHits, hitsHandle);
        //std::unordered_map< size_t, geo::WireID > hitToWire;
        //hitToWire.reserve(eventHits.size());
        
        const std::vector< art::Ptr< recob::Hit > > beamPFP_hits = pfpUtil.GetPFParticleHits_Ptrs( *particle, evt, fPFParticleTag );
        auto calo_dQdX = calo[index].dQdx();
        double reco_beam_len = thisTrack->Length();
        cout<<"$$Length"<<reco_beam_len<<endl;
        size_t i = 0;
        for (; thisTrack->Length(i) > reco_beam_len - sel_len; ++i ){
          ;
        }
        TVector3 pos = thisTrack->LocationAtPoint<TVector3>(i);
        cout<<"$$"<<"\tX"<<pos.X()<<"\tY"<<pos.Y()<<"\tZ"<<pos.Z()<<endl;
        
        auto TpIndices = calo[index].TpIndices();
        
        double wirecoord_U = fGeom->WireCoordinate(pos.Y(), pos.Z(), 0, 1, 0);
        double wirecoord_V = fGeom->WireCoordinate(pos.Y(), pos.Z(), 1, 1, 0);
        double wirecoord_X = fGeom->WireCoordinate(pos.Y(), pos.Z(), 2, 1, 0);
        cout<<"$$$WireCoord: U "<<wirecoord_U<<"\tV "<<wirecoord_V<<"\tX "<<wirecoord_X<<endl;
        
        std::vector< art::Ptr< recob::Hit > > remove_hits;
        //for (auto hit : beamPFP_hits){
        for (size_t kk = 0; kk < beamPFP_hits.size(); ++kk){
          auto hit = beamPFP_hits[kk];
          
          geo::WireID hitid = hit->WireID();
          //std::vector<geo::WireID> cwids = fGeom->ChannelToWire(hit->Channel());
          //cout<<"@@@"<<hitid.Plane<<"\t"<<hitid.Wire<<"\t"<<hitid.TPC<<"\t"<<hitid.Cryostat<<endl;
          
          if (hitid.TPC == 1 && hitid.Cryostat == 0) {
            if (hitid.Plane == 0) { // plane U
              if (hitid.Wire > wirecoord_U) {
                cout<<"$$hitid.Plane == 0"<<endl;
                //hcol.emplace_back(std::move(*hit), channelHitWires.at(kk));
                remove_hits.push_back(hit);
              }
            }
            else if (hitid.Plane == 1) { // plane V
              if (hitid.Wire < wirecoord_V) {
                cout<<"$$hitid.Plane == 1"<<endl;
                //hcol.emplace_back(std::move(*hit), channelHitWires.at(kk));
                remove_hits.push_back(hit);
              }
            }
            else if (hitid.Plane == 2) { // plane X
              if (hitid.Wire > wirecoord_X) {
                cout<<"$$hitid.Plane == 2"<<endl;
                //hcol.emplace_back(std::move(*hit), channelHitWires.at(kk));
                remove_hits.push_back(hit);
              }
            }
          }
        }
        // fill hits
        for (size_t kk = 0; kk < eventHits.size(); ++kk){
          auto hit = eventHits[kk];
          vector<art::Ptr<recob::Hit>>::iterator itr;
          itr = find(remove_hits.begin(), remove_hits.end(), hit);
          if (itr == remove_hits.end()) {
            hcol.emplace_back(std::move(*hit), channelHitWires.at(kk));
          }
        }
        
        cout<<"$$$reco_beam_len_sce "<<reco_beam_len_sce<<endl;
        
        
        /*auto TpIndices = calo[index].TpIndices();
        auto calo_dQdX = calo[index].dQdx();
        for( size_t i = 0; i < calo_dQdX.size(); ++i ){
          const recob::Hit & theHit = (*hitsHandle)[ TpIndices[i] ];
          cout<<"@@@TPC"<<theHit.WireID().TPC<<endl;
          cout<<"@@@Wire"<<theHit.WireID().Wire<<endl;
        }
        */
        
/*
        cout<<"@@@7"<<endl;
        auto assns = std::make_unique<art::Assns<recob::Hit, recob::SpacePoint>>();
        cout<<"@@@2"<<endl;
        std::vector< art::Ptr<recob::Hit> > eventHits;
        cout<<"@@@7"<<endl;
        art::fill_ptr_vector(eventHits, hitsHandle);
        cout<<"@@@2"<<endl;
        //art::FindManyP< recob::SpacePoint > spFromHit(hitsHandle, evt, fSpModuleLabel);
        cout<<"@@@7"<<endl;
        art::FindManyP< recob::Hit > hitsFromSp(spHandle, evt, fSpModuleLabel);
        cout<<"@@@2"<<endl;
        auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt);
        // map induction spacepoints to TPC by collection hits
        std::unordered_map< size_t, size_t > spToTPC;
        cout<<"@@@7"<<endl;
        for (size_t i = 0; i < spHandle->size(); ++i)
        {
            auto hits = hitsFromSp.at(i);
            size_t tpc = geo::WireID::InvalidID;
            for (const auto & h : hits) // find Collection hit, assume one is enough
            {
                if (h->SignalType() == geo::kCollection) { tpc = h->WireID().TPC; break; }
            }
            if (tpc == geo::WireID::InvalidID)
            {
        //mf::LogWarning("DisambigFromSpacePoints") << "No collection hit for this spacepoint.";
                continue;
            }
            for (const auto & h : hits) // set mapping for Induction hits
            {
                if (h->SignalType() == geo::kInduction) { spToTPC[i] = tpc; }
            }
        }
        cout<<"@@@3"<<endl;
        std::map< unsigned int, std::map< unsigned int, std::map< unsigned int, std::vector< size_t > > > > indHits;                        // induction hits resolved with spacepoints
        std::vector<size_t> unassignedHits;                   // hits to resolve by neighoring assignments
        std::unordered_map< size_t, geo::WireID > hitToWire;                 // final hit-wire assignments
        std::unordered_map< size_t, std::vector<geo::WireID> > hitToNWires;  // final hit-many-wires assignments

        hitToWire.reserve(eventHits.size());
*/
        /*int n = runOnSpacePoints(eventHits, spFromHit, spToTPC, hitToWire, indHits, unassignedHits);
        mf::LogInfo("DisambigFromSpacePoints") << n << " hits undisambiguated by space points.";

        if (fUseNeighbors)
        {
            n = resolveUnassigned(detProp, hitToWire, eventHits, indHits, unassignedHits, fNumNeighbors);
            mf::LogInfo("DisambigFromSpacePoints") << n << " hits undisambiguated by neighborhood.";
        }

        if (fMoveLeftovers == "repeat")     { assignEveryAllowedWire(hitToNWires, eventHits, unassignedHits);      }
        else if (fMoveLeftovers == "first") { assignFirstAllowedWire(hitToWire, eventHits, unassignedHits);        }
        else                { mf::LogInfo("DisambigFromSpacePoints") << "Remaining undisambiguated hits dropped."; }
        */
/*        auto const hitPtrMaker = art::PtrMaker<recob::Hit>(evt);
        cout<<"@@@4"<<endl;
        for (auto const & hw : hitToWire)
        {
            size_t key = hw.first;
            geo::WireID wid = hw.second;

            recob::HitCreator new_hit(*(eventHits[key]), wid);

            hcol.emplace_back(new_hit.move(), channelHitWires.at(key));//, channelHitRawDigits.at(key));

            auto hitPtr = hitPtrMaker(hcol.size() - 1);
            auto sps = spFromHit.at(eventHits[key].key());
            for (auto const & spPtr : sps)
            {
                assns->addSingle(hitPtr, spPtr);
            }
        }
        cout<<"@@@5"<<endl;
        for (auto const & hws : hitToNWires)
        {
            size_t key = hws.first;
            for (auto const & wid : hws.second)
            {
                recob::HitCreator new_hit(*(eventHits[key]), wid);

                hcol.emplace_back(new_hit.move(), channelHitWires.at(key));//, channelHitRawDigits.at(key));

                auto hitPtr = hitPtrMaker(hcol.size() - 1);
                auto sps = spFromHit.at(eventHits[key].key());
                for (auto const & spPtr : sps)
                {
                    assns->addSingle(hitPtr, spPtr);
                }
            }
        }
*/
        cout<<"@@@1"<<endl;
      }
    }
  }
  
  
  hcol.put_into(evt);
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
