////////////////////////////////////////////////////////////////////////
// Class:       MichelTiming
// Plugin Type: analyzer (art v3_05_01)
// File:        MichelTiming_module.cc
//
// Generated at Mon Jun 22 22:21:41 2020 by Tingjun Yang using cetskelgen
// from cetlib version v3_10_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RawData/RDTimeStamp.h"
#include "lardata/ArtDataHelper/MVAReader.h"

#include "TTree.h"
#include "TH1F.h"

#include <vector>

constexpr int kMaxWF = 2000;

using namespace std;

namespace pdsp {
  class MichelTiming;
}


class pdsp::MichelTiming : public art::EDAnalyzer {
public:
  explicit MichelTiming(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MichelTiming(MichelTiming const&) = delete;
  MichelTiming(MichelTiming&&) = delete;
  MichelTiming& operator=(MichelTiming const&) = delete;
  MichelTiming& operator=(MichelTiming&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;

private:

  //art::InputTag fOpDetModuleLabel;

  TTree *anatree;
  //TH1F *hist;
  int run;
  int event;
  vector<float> pandorat0;
  vector<int> trkid;
  vector<float> vtxx, vtxy, vtxz, endx, endy, endz;
  vector<float> michelscore;
  vector<int> michelhits;
  vector<short> pdchannel;
  vector<float> pdt0;
  vector<float> pd2t0;
  int nWF;
  int waveform[kMaxWF][2000];

  art::ServiceHandle<art::TFileService> tfs;

};


pdsp::MichelTiming::MichelTiming(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}//,
    //fOpDetModuleLabel(p.get< art::InputTag >("OpDetModuleLabel"))
{
}

void pdsp::MichelTiming::analyze(art::Event const& e)
{

  run = e.run();
  event = e.id().event();
  auto wfListHandle = e.getHandle< std::vector<raw::OpDetWaveform> >("ssprawdecoder:internal");
  auto wf2ListHandle = e.getHandle< std::vector<raw::OpDetWaveform> >("ssprawdecoder:external");

  art::FindOneP<raw::RDTimeStamp> assn(wfListHandle, e, "ssprawdecoder:internal");
  //cout<<assn.size()<<endl;
  art::FindOneP<raw::RDTimeStamp> assn2(wf2ListHandle, e, "ssprawdecoder:external");
  
  pdchannel.clear();
  pdt0.clear();
  nWF = 0;
  for (int i = 0; i<kMaxWF; ++i){
    for (int j = 0; j<2000; ++j){
      waveform[i][j] = 0;
    }
  }
  auto t0 = assn2.at(0)->GetTimeStamp();
  for (size_t i = 0; i<(*wfListHandle).size(); ++i){
    auto &wf = (*wfListHandle)[i];
//    auto &op1 = (*op1ListHandle)[i];
//    auto &op2 = (*op2ListHandle)[i];
    //for (auto const & wf : *wfListHandle){
    pdchannel.push_back(wf.ChannelNumber());
    auto t = assn.at(i)->GetTimeStamp();
    if (t>t0){
      pdt0.push_back((t-t0)/150.);
      //cout<<wf.ChannelNumber()<<" "<<t<<" "<<t0<<" "<<(t-t0)/150.<<endl;
    }
    else{
      pdt0.push_back(-((t0-t)/150.));
      //cout<<wf.ChannelNumber()<<" "<<t<<" "<<t0<<" "<<-((t0-t)/150.)<<endl;
    }
//    hist = tfs->make<TH1F>(Form("h%d_%d_%d", run, event, int(i)),
//                           Form("h%d_%d_%d", run, event, int(i)),
//                           2000,0,2000);
    if (nWF<kMaxWF){
      for(size_t j=0; j<wf.Waveform().size(); j++){
        waveform[nWF][j] =  wf.Waveform()[j];
      }
      ++nWF;
    }
  }
//  auto flashListHandle = e.getHandle < std::vector < recob::OpFlash > >(fFlashModuleLabel);
//  if (!flashListHandle) {
//    return;
//  }
//
//  std::size_t i = 0;
//  for (const recob::OpFlash & flash : *flashListHandle){
//    std::cout<< i++ <<" "<<flash.Time()<<" "<<flash.TotalPE()<<std::endl;
//  }

  pd2t0.clear();
  for (size_t i = 0; i<(*wf2ListHandle).size(); ++i){
    auto &wf = (*wf2ListHandle)[i];
//    auto &op1 = (*op1ListHandle)[i];
//    auto &op2 = (*op2ListHandle)[i];
    //for (auto const & wf : *wfListHandle){
    //pdchannel.push_back(wf.ChannelNumber());
    pd2t0.push_back(wf.TimeStamp());
  }
  
  pandorat0.clear();
  trkid.clear();
  vtxx.clear();
  vtxy.clear();
  vtxz.clear();
  endx.clear();
  endy.clear();
  endz.clear();
  michelscore.clear();
  michelhits.clear();
  auto t0ListHandle = e.getHandle< std::vector<anab::T0> >("pandora");
  if (!t0ListHandle) return;
  auto pfListHandle = e.getHandle< std::vector<recob::PFParticle> >("pandora");
  if (!pfListHandle) return;
  auto trkListHandle = e.getHandle< std::vector<recob::Track> >("pandoraTrack");
  if (!trkListHandle) return;
  std::vector < art::Ptr < recob::Track > > trkList;
  if (trkListHandle) {
    art::fill_ptr_vector(trkList, trkListHandle);
  }
  art::FindOneP<recob::PFParticle> assn3(t0ListHandle, e, "pandora");
  art::FindOneP<recob::Track> assn4(pfListHandle, e, "pandoraTrack");
  std::vector < art::Ptr < recob::Hit > > hitList;
  auto hitListHandle = e.getHandle < std::vector < recob::Hit > >("hitpdune");
  if (hitListHandle) {
    art::fill_ptr_vector(hitList, hitListHandle);
  }
  // Get track-hit association
  art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(trkListHandle, e,"pandoraTrack"); // to associate tracks and hits
  art::FindManyP<recob::Track> thass(hitListHandle, e, "pandoraTrack"); //to associate hit just trying
  anab::MVAReader<recob::Hit,4> hitResults(e, "emtrkmichelid:emtrkmichel");

  for (size_t i = 0; i<t0ListHandle->size(); ++i){
    auto &t0 = (*t0ListHandle)[i];
    if (assn3.at(i).isAvailable()){
      if (assn4.at(assn3.at(i).key()).isAvailable()){
        auto &trk = assn4.at(assn3.at(i).key());
        trkid.push_back(trk->ID());
        vtxx.push_back(trk->Vertex().X());
        vtxy.push_back(trk->Vertex().Y());
        vtxz.push_back(trk->Vertex().Z());
        endx.push_back(trk->End().X());
        endy.push_back(trk->End().Y());
        endz.push_back(trk->End().Z());
        pandorat0.push_back(t0.Time()*1e-3);
        int endwire = -1;
        int endtpc = -1;
        double endpeakt = -1;
        std::vector<int> wirekeys;
        if (fmthm.isValid()){
          float disend = 999999;
          auto vhit=fmthm.at(trk.key());
          auto vmeta=fmthm.data(trk.key());
          for (size_t ii = 0; ii<vhit.size(); ++ii){ //loop over all meta data hit
            bool fBadhit = false;
            if (vmeta[ii]->Index() == static_cast<unsigned int>(std::numeric_limits<int>::max())){
              fBadhit = true;
              continue;
            }
            if (vmeta[ii]->Index()>=trk->NumberTrajectoryPoints()){
              throw cet::exception("Calorimetry_module.cc") << "Requested track trajectory index "<<vmeta[ii]->Index()<<" exceeds the total number of trajectory points "<<trk->NumberTrajectoryPoints()<<" for track index "<<trkid.back()<<". Something is wrong with the track reconstruction. Please contact tjyang@fnal.gov!!";
            }
            if (!trk->HasValidPoint(vmeta[ii]->Index())){
              fBadhit = true;
              continue;
            }
            auto loc = trk->LocationAtPoint(vmeta[ii]->Index());
            if (fBadhit) continue; //HY::If BAD hit, skip this hit and go next
            if(vhit[ii]->WireID().Plane==2){
              wirekeys.push_back(vhit[ii].key());
              float dis = sqrt(pow(loc.X() - endx.back(), 2)+
                               pow(loc.Y() - endy.back(), 2)+
                               pow(loc.Z() - endz.back(), 2));
              if(dis<disend){
                disend = dis;
                endwire=vhit[ii]->WireID().Wire;
                endpeakt=vhit[ii]->PeakTime();
                endtpc=vhit[ii]->WireID().TPC;
              }
            }
          }
        }
        float avgmichelscore = 0;
        int nhits = 0;
        for(size_t hitl=0;hitl<hitList.size();hitl++){
          std::array<float,4> cnn_out=hitResults.getOutput(hitList[hitl]);
          auto & tracks = thass.at(hitList[hitl].key());
          // hit not on the track
          if (std::find(wirekeys.begin(), wirekeys.end(), hitList[hitl].key()) != wirekeys.end()) continue;
          // hit not on a long track
          if (!tracks.empty() && int(tracks[0].key()) != trkid.back() && trkList[tracks[0].key()]->Length()>25) continue;
          int planeid=hitList[hitl]->WireID().Plane;
          if (planeid!=2) continue;
          int tpcid=hitList[hitl]->WireID().TPC;
          if (tpcid!=endtpc) continue;
          float peakth1=hitList[hitl]->PeakTime();
          int wireh1=hitList[hitl]->WireID().Wire;
          if(std::abs(wireh1-endwire)<15 && std::abs(peakth1-endpeakt)<100 && tpcid==endtpc){
            ++nhits;
            avgmichelscore += cnn_out[hitResults.getIndex("michel")];
          }
        }
        if (nhits) avgmichelscore /= nhits;
        michelscore.push_back(avgmichelscore);
        michelhits.push_back(nhits);
      }
    }
  }
  
  anatree->Fill();

}

void pdsp::MichelTiming::beginJob()
{
  // Implementation of optional member function here.
  anatree = tfs->make<TTree>("anatree", "anatree");
  anatree->Branch("run", &run);
  anatree->Branch("event", &event);
  anatree->Branch("pandorat0", &pandorat0);
  anatree->Branch("trkid", &trkid);
  anatree->Branch("vtxx", &vtxx);
  anatree->Branch("vtxy", &vtxy);
  anatree->Branch("vtxz", &vtxz);
  anatree->Branch("endx", &endx);
  anatree->Branch("endy", &endy);
  anatree->Branch("endz", &endz);
  anatree->Branch("michelscore", &michelscore);
  anatree->Branch("michelhits", &michelhits);
  anatree->Branch("pdchannel", &pdchannel);
  anatree->Branch("pdt0", &pdt0);
  anatree->Branch("pd2t0", &pd2t0);
  anatree->Branch("nWF", &nWF, "nWF/I");
  anatree->Branch("waveform", waveform, "waveform[nWF][2000]/I");
}

DEFINE_ART_MODULE(pdsp::MichelTiming)
