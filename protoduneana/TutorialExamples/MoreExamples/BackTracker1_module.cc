////////////////////////////////////////////////////////////////////////
// Class:       BackTracker1
// Plugin Type: analyzer (Unknown Unknown)
// File:        BackTracker1_module.cc
//
// Generated at Sun Apr 10 22:23:05 2022 by Tingjun Yang using cetskelgen
// from  version .
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
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "TH1D.h"

using namespace std;

namespace pdsp {
  class BackTracker1;
}


class pdsp::BackTracker1 : public art::EDAnalyzer {
public:
  explicit BackTracker1(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  BackTracker1(BackTracker1 const&) = delete;
  BackTracker1(BackTracker1&&) = delete;
  BackTracker1& operator=(BackTracker1 const&) = delete;
  BackTracker1& operator=(BackTracker1&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  void beginJob() override;

private:

  // distance between reconstructed trajectory points and true energy deposition location
  TH1D *hdeltax;
  TH1D *hdeltay;
  TH1D *hdeltaz;

  // number of points in 14 cm segments
  TH1D *hsegpts_reco; //using reconstructed information
  TH1D *hsegpts_true; //using simchannel information
  TH1D *hsegpts_sed; //using simenergydeposit information

};


pdsp::BackTracker1::BackTracker1(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void pdsp::BackTracker1::analyze(art::Event const& e)
{
  // Get services
  auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;

  // get tracks
  std::vector<art::Ptr<recob::Track> > tracklist;
  auto trackListHandle = e.getHandle< std::vector<recob::Track> >("pandoraTrack");
  if (trackListHandle){
    art::fill_ptr_vector(tracklist, trackListHandle);
  }
  else{
    cout<<"trackListHandle invalid"<<endl;
    return;
  }

  art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(trackListHandle, e, "pandoraTrack");
  if (!fmthm.isValid()){
    cout<<"fmthm invalid"<<endl;
  }

  for (size_t trkIter = 0; trkIter < tracklist.size(); ++trkIter) {

    auto & track = tracklist[trkIter];

    //Ignore short tracks
    if (track->Length()<10) continue;
    
    auto vhit = fmthm.at(trkIter);
    auto vmeta = fmthm.data(trkIter);

    int nsegpts_reco = 0;
    int nsegpts_true = 0;
    double seglength_reco = 0;
    double seglength_true = 0;
    double prevx_reco = 0;
    double prevy_reco = 0;
    double prevz_reco = 0;
    double prevx_true = 0;
    double prevy_true = 0;
    double prevz_true = 0;
    bool isfirst_reco = true;
    bool isfirst_true = true;
    for (size_t ip = 0; ip<track->NPoints(); ++ip){
      // Ignore invalid trajectory points
      if (!track->HasValidPoint(ip)) continue;
      auto &loc = track->LocationAtPoint(ip);
      //find number of points in reco segments of 14 cm
      if (isfirst_reco){
        prevx_reco = loc.X();
        prevy_reco = loc.Y();
        prevz_reco = loc.Z();
        nsegpts_reco = 1;
        isfirst_reco = false;
      }
      else{
        double dis = sqrt(pow(loc.X() - prevx_reco, 2)+
                          pow(loc.Y() - prevy_reco, 2)+
                          pow(loc.Z() - prevz_reco, 2));
        seglength_reco += dis;
        if (seglength_reco>14){
          hsegpts_reco->Fill(nsegpts_reco);
          //cout<<"Fill reco segment "<<nsegpts_reco<<endl;
          nsegpts_reco = 0;
          seglength_reco = 0;
        }
        ++nsegpts_reco;
        prevx_reco = loc.X();
        prevy_reco = loc.Y();
        prevz_reco = loc.Z();
      }
      //cout<<"reco x "<<loc.X()<<" reco y "<<loc.Y()<<" reco z "<<loc.Z()<<" seglength_reco "<<seglength_reco<<" nsegpts_reco "<<nsegpts_reco<<endl;
      for (size_t ii = 0; ii < vhit.size(); ++ii) {
        // Look for hit corresponding to the trajectory point
        if (vmeta[ii]->Index() == ip){
          auto simIDEs = bt_serv->HitToSimIDEs_Ps(clock_data, vhit[ii]);
          if (!simIDEs.empty()){
            // Find true location of the hit
//            for (auto const& ide : simIDEs){
//              cout<<"id = "<<ide->trackID<<" "<<ide->numElectrons<<" "<<ide->x<<" "<<ide->y<<" "<<ide->z<<" "<<endl;
//            }
            auto xyz = bt_serv->HitToXYZ(clock_data, vhit[ii]);
            //cout<<loc.X()<<" "<<loc.Y()<<" "<<loc.Z()<<" "<<xyz[0]<<" "<<xyz[1]<<" "<<xyz[2]<<" "<<string(vhit[ii]->WireID())<<endl;
            hdeltax->Fill(loc.X() - xyz[0]);
            hdeltay->Fill(loc.Y() - xyz[1]);
            hdeltaz->Fill(loc.Z() - xyz[2]);
            //find number of points in reco segments of 14 cm
            if (isfirst_true){
              prevx_true = xyz[0];
              prevy_true = xyz[1];
              prevz_true = xyz[2];
              nsegpts_true = 1;
              isfirst_true = false;
            }
            else{
              double dis = sqrt(pow(xyz[0] - prevx_true, 2)+
                                pow(xyz[1] - prevy_true, 2)+
                                pow(xyz[2] - prevz_true, 2));
              seglength_true += dis;
              if (seglength_true>14){
                hsegpts_true->Fill(nsegpts_true);
                //cout<<"Fill true segment "<<nsegpts_true<<endl;
                nsegpts_true = 0;
                seglength_true = 0;
              }
              ++nsegpts_true;
              prevx_true = xyz[0];
              prevy_true = xyz[1];
              prevz_true = xyz[2];
            }
            //cout<<"true x "<<xyz[0]<<" true y "<<xyz[1]<<" true z "<<xyz[2]<<" seglength_true "<<seglength_true<<" nsegpts_true "<<nsegpts_true<<endl;
          }
        }
      }
    }
  }

  // Get SimEnergyDeposit
  std::vector<art::Ptr<sim::SimEnergyDeposit>> sedlist;
  auto sedListHandle = e.getHandle< std::vector<sim::SimEnergyDeposit> >("largeant:LArG4DetectorServicevolTPCActive");
  if (sedListHandle){
    art::fill_ptr_vector(sedlist, sedListHandle);
  }
  else{
    cout<<"sedListHandle invalid"<<endl;
    return;
  }

  int nsegpts_sed = 0;
  double seglength_sed = 0;
  double prevx_sed = 0;
  double prevy_sed = 0;
  double prevz_sed = 0;
  bool isfirst_sed = true;
  for (size_t i = 0; i<sedlist.size(); ++i){
    double x = sedlist[i]->X();
    double y = sedlist[i]->Y();
    double z = sedlist[i]->Z();
    int trkid = sedlist[i]->TrackID();
    //int pdg = sedlist[i]->PdgCode();
    //cout<<x<<" "<<y<<" "<<z<<" "<<trkid<<" "<<pdg<<endl;
    if (trkid!=1) continue;
    if (isfirst_sed){
      prevx_sed = x;
      prevy_sed = y;
      prevz_sed = z;
      nsegpts_sed = 1;
      isfirst_sed = false;
    }
    else{
      double dis = sqrt(pow(x - prevx_sed, 2)+
                        pow(y - prevy_sed, 2)+
                        pow(z - prevz_sed, 2));
      seglength_sed += dis;
      if (seglength_sed>14){
        hsegpts_sed->Fill(nsegpts_sed);
        //cout<<"Fill sed segment "<<nsegpts_sed<<endl;
//        if (nsegpts_sed<100){
//          cout<<x<<" "<<y<<" "<<z<<" "<<prevx_sed<<" "<<prevy_sed<<" "<<prevz_sed<<" "<<nsegpts_sed<<endl;
//        }
        nsegpts_sed = 0;
        seglength_sed = 0;
      }
      ++nsegpts_sed;
      prevx_sed = x;
      prevy_sed = y;
      prevz_sed = z;
    }
    //cout<<"sed x "<<x<<" sed y "<<y<<" sed z "<<z<<" seglength_sed "<<seglength_sed<<" nsegpts_sed "<<nsegpts_sed<<endl;
  }
}

void pdsp::BackTracker1::beginJob(){

  art::ServiceHandle< art::TFileService > tfs;
  hdeltax = tfs->make<TH1D>("hdeltax", ";Reco x - true x (cm);N trajectory points",100,-2,2);
  hdeltay = tfs->make<TH1D>("hdeltay", ";Reco y - true y (cm);N trajectory points",100,-2,2);
  hdeltaz = tfs->make<TH1D>("hdeltaz", ";Reco z - true z (cm);N trajectory points",100,-2,2);
  hsegpts_reco = tfs->make<TH1D>("hsegpts_reco", "Reco segments;Number of points;Segments", 200,0,200);
  hsegpts_true = tfs->make<TH1D>("hsegpts_true", "True segments;Number of points;Segments", 200,0,200);
  hsegpts_sed = tfs->make<TH1D>("hsegpts_sed", "SimEnergyDeposit segments;Number of points;Segments", 1000,0,1000);
}

DEFINE_ART_MODULE(pdsp::BackTracker1)
