////////////////////////////////////////////////////////////////////////
// Class:       CheckHitsAndTracks
// Plugin Type: analyzer (art v3_05_01)
// File:        CheckHitsAndTracks_module.cc
//
// Generate simple summary tree from tracks and hits to check sim and reco
// Follow the same logic as in CosmicsdQdx
//
// Generated at 31/01/2023 by Thibaut Houdy
//////////////////////////////////////////////////////////////////////Check//
#include <iostream>
#include <vector>
#include <utility>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "art_root_io/TFileService.h"

// LArSoft includes
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/Exceptions.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/ArtDataHelper/TrackUtils.h"

#include "TPolyMarker3D.h"
#include "TPolyMarker.h"
#include "TPolyLine.h"
#include "TPolyLine3D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TLine.h"


#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"

// ROOT
#include "TTree.h"

//
using std::cerr;
using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::pair;

//
namespace pdvdana {
  class CheckHitsAndTracks;
}

namespace {

  class dist_projected{
  public:
    dist_projected(recob::Hit const& h, geo::Geometry const* g):
      hit(h), geom(g){}
    bool operator() (std::pair<geo::WireID,float> const& i, std::pair<geo::WireID,float> const& j)
    {
      float dw_i = ((int)(i.first.Wire) - (int)(hit.WireID().Wire))*geom->WirePitch(i.first.asPlaneID());
      float dw_j = ((int)(j.first.Wire) - (int)(hit.WireID().Wire))*geom->WirePitch(j.first.asPlaneID());
      float dt_i = i.second - hit.PeakTime();
      float dt_j = j.second - hit.PeakTime();
      return (std::sqrt(dw_i*dw_i + dt_i*dt_i) < std::sqrt(dw_j*dw_j + dt_j*dt_j));
    }
  private:
    recob::Hit const& hit;
    geo::Geometry const* geom;
  };
}


class pdvdana::CheckHitsAndTracks : public art::EDAnalyzer {
public:
  explicit CheckHitsAndTracks(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CheckHitsAndTracks(CheckHitsAndTracks const&) = delete;
  CheckHitsAndTracks(CheckHitsAndTracks&&) = delete;
  CheckHitsAndTracks& operator=(CheckHitsAndTracks const&) = delete;
  CheckHitsAndTracks& operator=(CheckHitsAndTracks&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  int      fLogLevel;
  int      fFlagInfos = 1;
  int      fFlagWarning = 2;
  int      fFlagDebug = 3;

  string   fHitModuleLabel;
  string   fTrackModuleLabel;
  float    fTrackMinLen;
  int      ev_num;

  unsigned fTotalTracks;
  unsigned fSelectedTracks;

  // summary tree
  TTree   *ftrackTree;
  TTree   *frecoTree;

  TCanvas* canvas_xy; 
  TCanvas* canvas_3D;
  TCanvas* canvas_tpc;
  
  //
  unsigned fEventNum;
  unsigned fTrackId;
  unsigned fTrajPoints;
  unsigned fNtpcs, fNplanes;
  unsigned ftrackStartTick;
  unsigned ftrackEndTick;
  float ftrackLen;
  float ftrackDx;
  float ftrackDy;
  float ftrackDz;
  float ftrackStartX;
  float ftrackStartY;
  float ftrackStartZ;
  float ftrackEndX;
  float ftrackEndY;
  float ftrackEndZ;
  float ftracktheta;
  float ftrackphi;
  float ftracknorm;

  float frecoX;
  float frecoY;
  float frecoZ;
  unsigned frecoPlane, frecoTPC, frecoWire, frecoChannel;
  int frecoTime;

  vector<int> fOffsetWireID_u;
  vector<int> fOffsetWireID_v;
  vector<int> fOffsetWireID_z;
  
  float    fgeoXmin = 1e6; 
  float    fgeoXmax =-1e6; 
  float    fgeoYmin = 1e6; 
  float    fgeoYmax =-1e6; 
  float    fgeoZmin = 1e6; 
  float    fgeoZmax =-1e6; 

  // detector geometry
  const geo::Geometry* fGeom;

  void checkTrackCharacs( float Dx, float Dy, float Dz, vector<float> &track_char );
  void DrawSquare(TPolyLine *Boite, float a_1, float a_2, float b_1, float b_2);  
  void DrawCube(TPolyLine3D *Boite, float x_min, float y_min, float z_min, float x_max, float y_max, float z_max );
  void SetHitStyle(TPolyMarker *l_hits);
  unsigned GetWireOffset(unsigned plane_id, unsigned tpc_id);

};


pdvdana::CheckHitsAndTracks::CheckHitsAndTracks(fhicl::ParameterSet const& p)
  : EDAnalyzer{p} ,
  fLogLevel( p.get< int >("LogLevel") ),
  fHitModuleLabel( p.get< std::string >("HitModuleLabel") ),
  fTrackModuleLabel( p.get< std::string >("TrackModuleLabel") ),
  fTrackMinLen( p.get< float  >("TrackMinLen") )
  {
    fGeom    = &*art::ServiceHandle<geo::Geometry>();
  }

//
void pdvdana::CheckHitsAndTracks::analyze(art::Event const& ev)
{

  const string myname = "pdvdana::CheckHitsAndTracks::analyze: ";

  // get tracks

  string hittag = "gaushit";
  art::InputTag hit_tag(fHitModuleLabel);
  string trcktag = "pandoraTrack";
  art::InputTag trck_tag(trcktag);
  string hittag_dis = "hitpdune";
  art::InputTag hit_tag_dis(hittag_dis);
  string calotag = "pandoraGnocalo";//"calo";
  art::InputTag calo_tag(calotag);
  string clutag = "pandora";
  art::InputTag cluster_tag(clutag);

  vector<vector<float>> trackSpacePoints;
  vector<vector<float>> vertexSpacePoints;
  vector<vector<int>> vector_Hits_tpc_channel(fNplanes);
  vector<vector<int>> vector_Hits_tpc_time(fNplanes);

  vector<vector<int> > fTrackStartWireID(fNplanes);
  vector<vector<int> > fTrackEndWireID(fNplanes);

  vector<TPolyMarker*> EndTracks_tpc(fNplanes);
  vector<TPolyMarker*> StartTracks_tpc(fNplanes);
  vector<TPolyMarker*> Hits_tpc(fNplanes);

  ev_num++;

  auto Tracks = ev.getValidHandle<vector<recob::Track>>(trck_tag);
  auto Hits = ev.getValidHandle<vector<recob::Hit>>(hit_tag);

  const art::FindManyP<recob::Track> fmp(Hits ,ev ,trcktag);

  fTotalTracks += Tracks->size();
  fEventNum = ev_num;

  vector<float> track_charac(6);

  auto const hitHandle = ev.getValidHandle<std::vector<recob::Hit>>(hit_tag);
  std::vector<art::Ptr<recob::Hit>> hits;
  art::fill_ptr_vector(hits, hitHandle);
  
  if (fLogLevel>=fFlagInfos){
    cout<<myname<<Hits->size()<<" hits found"<<endl;
    cout<<myname<<Tracks->size()<<" tracks found"<<endl;
  }

  // Get the hits associated with the space points
  const art::FindManyP<recob::SpacePoint> fmsph(hitHandle, ev, cluster_tag);
  if (!fmsph.isValid()) {
    throw cet::exception("LArPandoraShowerCheatingAlg")
      << "Spacepoint and hit association not valid. Stopping.";
  }

  //----------------------------------------------------------------------------------------------
  //Loops on tracks
  for (unsigned itrk = 0; itrk < Tracks->size(); ++itrk) {
    const recob::Track& track = Tracks->at(itrk);
    ftrackLen    = track.Length();     

    if (ftrackLen<fTrackMinLen)
      continue;
    ftrackStartX  = track.Start().X();
    ftrackStartY  = track.Start().Y();
    ftrackStartZ  = track.Start().Z();
    ftrackStartTick  = 0;    
    ftrackEndX    = track.End().X();
    ftrackEndY    = track.End().Y();
    ftrackEndZ    = track.End().Z();
    ftrackEndTick = 0;

    // Filling the track tree
    checkTrackCharacs(ftrackEndX-ftrackStartX, ftrackEndY-ftrackStartY, ftrackEndZ-ftrackStartZ, track_charac);
    ftrackDx = track_charac[0];
    ftrackDy = track_charac[1];
    ftrackDz = track_charac[2];
    ftracknorm = track_charac[3];
    ftracktheta = track_charac[4];
    ftrackphi = track_charac[5];
    fTrackId     = 0;
    ftrackTree->Fill();

    if (fLogLevel>=fFlagInfos){
      cout <<myname<< "-- Track: "<<itrk<<" x,y,z: (" <<ftrackStartX<<"; "<< ftrackStartY<<"; "<< ftrackStartZ<<" ) | x,y,z: ("<< ftrackEndX<<"; "<< ftrackEndY<<"; "<< ftrackEndZ<<")"<<endl;
    }

    geo::TPCID tpc_s = fGeom->FindTPCAtPosition(track.Start());
    geo::TPCID tpc_e = fGeom->FindTPCAtPosition(track.End());

    int wireID_s = -9999;
    int wireID_e = -9999;

    if(! (tpc_s.isValid && tpc_e.isValid)){
      if (fLogLevel>=fFlagWarning){
        cout<<"---!!--- TPCs are not valid : "<<endl;
        cout<<"---!!--- TPCs ID: "<<tpc_s.TPC <<" and "<<tpc_e.TPC << endl;
      } 
      continue;
    }
    else{
      geo::WireID g_wireID_s;
      geo::WireID g_wireID_e;
      for(unsigned int id_pl = 0; id_pl<fNplanes; id_pl++){
        if ((tpc_s.TPC != tpc_e.TPC) && (fLogLevel>=fFlagWarning)){
          cout<<"---> Carefull, start and end are in different TPCs"<<endl;
        }
        try{
          g_wireID_s = fGeom->NearestWireID(track.Start(), geo::PlaneID(0, tpc_s.TPC, id_pl));
        }
        catch(geo::InvalidWireError const& e_1) {
          g_wireID_s = e_1.suggestedWireID(); // pick the closest valid wire
        }
        //Looking for the wire associated with the end of the track
        try{
          g_wireID_e = fGeom->NearestWireID(track.End(), geo::PlaneID(0, tpc_e.TPC, id_pl));
        }
        catch(geo::InvalidWireError const& e_2) {
          g_wireID_e = e_2.suggestedWireID(); // pick the closest valid wire
        }

        wireID_s = g_wireID_s.Wire;
        wireID_e = g_wireID_e.Wire;
  
        if ((wireID_e!=-9999) && (wireID_s!=-9999)){

          wireID_s += GetWireOffset(id_pl, tpc_s.TPC);
          wireID_e += GetWireOffset(id_pl, tpc_e.TPC);
          if (fLogLevel>=fFlagDebug)
            cout<<"---- START - TPC: "<<tpc_s.TPC <<" wire ID: "<< wireID_s<<" | END : "<<tpc_e.TPC << " wire ID: "<< wireID_e << endl;        
          fTrackStartWireID[id_pl].push_back(g_wireID_s.Wire + GetWireOffset(id_pl, tpc_s.TPC));
          fTrackEndWireID[id_pl].push_back( g_wireID_e.Wire + GetWireOffset(id_pl, tpc_e.TPC)); 
          fSelectedTracks++;
        }    
      }
    }
    vector<float> pos_track = {ftrackStartX, ftrackStartY, ftrackStartZ, ftrackEndX, ftrackEndY, ftrackEndZ};
    trackSpacePoints.push_back(pos_track);
  }// Tracks loop

  //----------------------------------------------------------------------------------------------
  //Loop on hits
  for (auto hit : hits) {
    std::vector<art::Ptr<recob::Track>> tracks = fmp.at(hit.key());
    std::vector<art::Ptr<recob::SpacePoint>> sps = fmsph.at(hit.key());

    if (sps.size() == 1 && tracks.size()==1) {
      art::Ptr<recob::SpacePoint> sp = sps.front();
      frecoX = sp->XYZ()[0];
      frecoY = sp->XYZ()[1];
      frecoZ = sp->XYZ()[2];
      fTrackId     = 0;

      geo::WireID hit_wireID = hit->WireID();
      frecoTPC = hit_wireID.TPC;
      frecoPlane = hit_wireID.Plane;
      frecoWire = hit_wireID.Wire+GetWireOffset(frecoPlane, frecoTPC);
      frecoTime = hit->PeakTime();
      frecoChannel = hit->Channel();
      if (frecoTime>10000)
        continue;

      std::vector<geo::WireID> cwids = fGeom->ChannelToWire(frecoChannel);
      vector_Hits_tpc_channel[frecoPlane].push_back(frecoWire);
      vector_Hits_tpc_time[frecoPlane].push_back(frecoTime);      
      vector<float> pos_sp = {frecoX, frecoY, frecoZ};
      vertexSpacePoints.push_back(pos_sp);
      frecoTree->Fill();
      } 
  }// Hits loop

  //----------------------------------------------------------------------------------------------
  // Drawing
  for(unsigned int id_pl =0; id_pl <fNplanes; id_pl++){
    int max_channel=0;
    int max_tick = 0;

    StartTracks_tpc[id_pl] = new TPolyMarker(fTrackStartWireID[id_pl].size());
    EndTracks_tpc[id_pl] = new TPolyMarker(fTrackEndWireID[id_pl].size());
    Hits_tpc[id_pl] = new TPolyMarker(vector_Hits_tpc_channel[id_pl].size());
    //Setting the track polymarker
    for(unsigned int j=0; j<fTrackStartWireID[id_pl].size(); j++){
      StartTracks_tpc[id_pl]->SetPoint(j, fTrackStartWireID[id_pl][j], 2000);
      EndTracks_tpc[id_pl]->SetPoint(j, fTrackEndWireID[id_pl][j], 4000);
    }
    //Setting the hits polymarker
    for(unsigned int j=0; j<vector_Hits_tpc_channel[id_pl].size(); j++){
      Hits_tpc[id_pl]->SetPoint(j, vector_Hits_tpc_channel[id_pl][j], vector_Hits_tpc_time[id_pl][j]);
      if (fLogLevel >= fFlagDebug)
        cout<<"hits : "<<id_pl<<"  "<<j<<"  "<< vector_Hits_tpc_channel[id_pl][j]<<"  "<<vector_Hits_tpc_time[id_pl][j]<<endl;
      if (max_tick<vector_Hits_tpc_time[id_pl][j])  max_tick = vector_Hits_tpc_time[id_pl][j];
      if (max_channel<vector_Hits_tpc_channel[id_pl][j])  max_channel = vector_Hits_tpc_channel[id_pl][j];         
    }
    //cout<<"--------------- Plane : "<<id_pl<<"  "<<max_tick<<"  "<<max_channel<<endl;
  }

  TPolyMarker3D* otherPoly = new TPolyMarker3D(vertexSpacePoints.size());
  TPolyMarker *HitsYZ = new TPolyMarker(vertexSpacePoints.size());
  TPolyMarker *HitsXY = new TPolyMarker(vertexSpacePoints.size());
  TPolyMarker *HitsXZ = new TPolyMarker(vertexSpacePoints.size());

  int point=0;

  for (auto const& sp : vertexSpacePoints) {
    otherPoly->SetPoint(point, sp[0], sp[1], sp[2]);
    HitsYZ->SetPoint(point, sp[1], sp[2]);
    HitsXY->SetPoint(point, sp[0], sp[1]);
    HitsXZ->SetPoint(point, sp[0], sp[2]);
    ++point;
  }

  SetHitStyle(HitsYZ);
  SetHitStyle(HitsXY);
  SetHitStyle(HitsXZ);
  otherPoly->SetMarkerStyle(8);
  otherPoly->SetMarkerSize(1);    
  otherPoly->SetMarkerColor(kBlue);

  TPolyLine3D* Boite3D = new TPolyLine3D(17);
  TPolyLine *BoiteYZ = new TPolyLine(5);
  TPolyLine *BoiteXY = new TPolyLine(5);
  TPolyLine *BoiteXZ = new TPolyLine(5);
  DrawSquare(BoiteXZ, fgeoXmin, fgeoXmax, fgeoZmin, fgeoZmax);
  DrawSquare(BoiteYZ, fgeoXmin, fgeoYmax, fgeoZmin, fgeoZmax);
  DrawSquare(BoiteXY, fgeoXmin, fgeoXmax, fgeoYmin, fgeoYmax);
  DrawCube(Boite3D, fgeoXmin, fgeoYmin, fgeoZmin, fgeoXmax, fgeoYmax, fgeoZmax);

  art::ServiceHandle<art::TFileService> tfs;

  //Drawing
  TString canvasName_3D = Form("canvas3D_%i", ev_num);
  canvas_3D = tfs->make<TCanvas>(canvasName_3D, canvasName_3D);

  gStyle->SetOptStat(0);
  canvas_3D->cd();
  TH3F *axes = new TH3F("axes", "3D Hits and tracks distribution; X; Y; Z", 1, -400, 400, 1, -400, 400, 1, -400, 400);
  axes->SetDirectory(0);
  axes->Draw();
  Boite3D->Draw();
  otherPoly->Draw();   

  for (auto const& sp : trackSpacePoints) {
    TPolyLine3D* t_line = new TPolyLine3D(2);
    t_line->SetPoint(0, sp[0], sp[1], sp[2]);
    t_line->SetPoint(1, sp[3], sp[4], sp[5]);
    t_line->SetLineColor(2);
    t_line->SetLineWidth(2);
    t_line->Draw();
  }

  TString canvasName_xy = Form("canvasxy_%i", ev_num);
  canvas_xy = tfs->make<TCanvas>(canvasName_xy, canvasName_xy);
  TH2F *axos_XY = new TH2F("axes_1", "Y(X); X; Y", 1, -400, 400, 1, -400, 400);
  TH2F *axos_XZ = new TH2F("axes_2", "Z(X); X; Z", 1, -400, 400, 1, -400, 400);
  TH2F *axos_YZ = new TH2F("axes_3", "Z(Y); Y; Z", 1, -400, 400, 1, -400, 400);
  canvas_xy->Divide(2,2);
  canvas_xy->cd(1);
  axos_XY->Draw();
  BoiteXY->Draw();
  HitsXY->Draw();
  canvas_xy->cd(2);
  axos_XZ->Draw();    
  BoiteXZ->Draw();
  HitsXZ->Draw();
  canvas_xy->cd(3);
  axos_YZ->Draw();
  BoiteYZ->Draw();
  HitsYZ->Draw();

  TString vertexName = Form("space_points_%i", ev_num);
  for (auto const& sp : trackSpacePoints) {
    TLine* l1 = new TLine(sp[0], sp[1], sp[3], sp[4]);
    TLine* l2 = new TLine(sp[0], sp[2], sp[3], sp[5]);
    TLine* l3 = new TLine(sp[1], sp[2], sp[4], sp[5]);
    canvas_xy->cd(1);
    l1->SetLineColor(2);
    l1->SetLineWidth(1);
    l1->Draw();
    canvas_xy->cd(2);
    l2->SetLineColor(2);
    l2->SetLineWidth(1);
    l2->Draw();
    canvas_xy->cd(3);
    l3->SetLineColor(2);
    l3->SetLineWidth(1);
    l3->Draw();
  }

  TString canvasName_tpc = Form("canvastpc_%i", ev_num);
  canvas_tpc= tfs->make<TCanvas>(canvasName_tpc, canvasName_tpc);
  TH2F *axos_tpcU = new TH2F("U plane", "U plane; wire; time [ticks]", 1, 0, 4600, 1, 0, 6000);
  TH2F *axos_tpcV = new TH2F("V plane", "V plane; wire; time [ticks]", 1, 0, 4600, 1, 0, 6000);
  TH2F *axos_tpcZ = new TH2F("Z plane", "Z plane; wire; time [ticks]", 1, 0, 4600, 1, 0, 6000);

  for (unsigned i=0; i<fNplanes; i++ ){
    SetHitStyle(Hits_tpc[i]);
    StartTracks_tpc[i]->SetMarkerStyle(23);
    StartTracks_tpc[i]->SetMarkerSize(1);    
    StartTracks_tpc[i]->SetMarkerColor(kRed);
    EndTracks_tpc[i]->SetMarkerStyle(22);
    EndTracks_tpc[i]->SetMarkerSize(1);    
    EndTracks_tpc[i]->SetMarkerColor(kBlack);
  }
  SetHitStyle(Hits_tpc[1]);
  SetHitStyle(Hits_tpc[2]);
  
  canvas_tpc->Divide(2,2);
  canvas_tpc->cd(1);
  axos_tpcU->Draw();  
  EndTracks_tpc[0]->Draw();
  StartTracks_tpc[0]->Draw();
  Hits_tpc[0]->Draw();
  for (unsigned i =0; i < fNtpcs; i++){
    TLine* lines_U = new TLine(fOffsetWireID_u[i],0,fOffsetWireID_u[i],6000);
    canvas_tpc->cd(1);
    lines_U->SetLineColor(12);
    lines_U->SetLineWidth(1);
    lines_U->Draw();

  }

  canvas_tpc->cd(2);
  axos_tpcV->Draw();    
  EndTracks_tpc[1]->Draw();
  StartTracks_tpc[1]->Draw();
  Hits_tpc[1]->Draw();
  for (unsigned i =0; i < fNtpcs; i++){
    TLine* lines_V = new TLine(fOffsetWireID_v[i],0,fOffsetWireID_v[i],6000);
    canvas_tpc->cd(2);
    lines_V->SetLineColor(12);
    lines_V->SetLineWidth(1);
    lines_V->Draw();
  }

  canvas_tpc->cd(3);
  axos_tpcZ->Draw();
  EndTracks_tpc[2]->Draw();
  StartTracks_tpc[2]->Draw();
  Hits_tpc[2]->Draw();
  for (unsigned i =0; i < fNtpcs; i++){
    TLine* lines_Z = new TLine(fOffsetWireID_z[i],0,fOffsetWireID_z[i],6000);
    canvas_tpc->cd(3);
    lines_Z->SetLineColor(12);
    lines_Z->SetLineWidth(1);
    lines_Z->Draw();
  }

  canvas_tpc->Write(canvasName_tpc);
  canvas_3D->Write(canvasName_3D);
  canvas_xy->Write(canvasName_xy);


} 
//
void pdvdana::CheckHitsAndTracks::beginJob()
{

  ev_num = 0;
  fTotalTracks    = 0;
  fSelectedTracks = 0;
  fNplanes = fGeom->Nplanes();
  fNtpcs = fGeom->NTPC();

  fOffsetWireID_u.clear();
  fOffsetWireID_v.clear();
  fOffsetWireID_z.clear();

  unsigned int mem_nwires_u = 0;
  unsigned int mem_nwires_v = 0;
  unsigned int mem_nwires_z = 0;
  cout<<"Geometry : "<<endl;
  cout<<"number of planes :  "<<fNplanes<<endl;
  cout<<"number of tpc :  "<<fNtpcs<<endl;

  for(unsigned t_tpc_id=0;t_tpc_id<fNtpcs;t_tpc_id++){

    geo::TPCID tpcid{{0}, t_tpc_id};
    geo::PlaneID const uplane_id{tpcid, geo::View_t::kU};
    geo::PlaneID const vplane_id{tpcid, geo::View_t::kV};
    geo::PlaneID const zplane_id{tpcid, geo::View_t::kZ};
    unsigned int nwires_u = fGeom->Nwires(uplane_id);
    unsigned int nwires_v = fGeom->Nwires(vplane_id);
    unsigned int nwires_z = fGeom->Nwires(zplane_id);

    fOffsetWireID_u.push_back(mem_nwires_u);
    fOffsetWireID_v.push_back(mem_nwires_v);
    fOffsetWireID_z.push_back(mem_nwires_z);

    mem_nwires_u += nwires_u;
    mem_nwires_v += nwires_v;
    mem_nwires_z += nwires_z;


    std::cout << "  TPC " << t_tpc_id << " center: ("<< fGeom->TPC(tpcid).GetCenter().X()      << "," << fGeom->TPC(tpcid).GetCenter().Y()      << ","<< fGeom->TPC(tpcid).GetCenter().Z() << ")"
                               << " box:  ["  << fGeom->TPC(tpcid).BoundingBox().MinX() << "," << fGeom->TPC(tpcid).BoundingBox().MaxX() << "]" 
                                      << "["  << fGeom->TPC(tpcid).BoundingBox().MinY() << "," << fGeom->TPC(tpcid).BoundingBox().MaxY() << "]"
                                      << "["  << fGeom->TPC(tpcid).BoundingBox().MinZ() << "," << fGeom->TPC(tpcid).BoundingBox().MaxZ() << "]" << std::endl;

    if(fgeoXmin > fGeom->TPC(tpcid).BoundingBox().MinX()) fgeoXmin = fGeom->TPC(tpcid).BoundingBox().MinX();
    if(fgeoXmax < fGeom->TPC(tpcid).BoundingBox().MaxX()) fgeoXmax = fGeom->TPC(tpcid).BoundingBox().MaxX();
    if(fgeoYmin > fGeom->TPC(tpcid).BoundingBox().MinY()) fgeoYmin = fGeom->TPC(tpcid).BoundingBox().MinY();
    if(fgeoYmax < fGeom->TPC(tpcid).BoundingBox().MaxY()) fgeoYmax = fGeom->TPC(tpcid).BoundingBox().MaxY();
    if(fgeoZmin > fGeom->TPC(tpcid).BoundingBox().MinZ()) fgeoZmin = fGeom->TPC(tpcid).BoundingBox().MinZ();
    if(fgeoZmax < fGeom->TPC(tpcid).BoundingBox().MaxZ()) fgeoZmax = fGeom->TPC(tpcid).BoundingBox().MaxZ();
    cout<< "Number of wires : "<<nwires_u<<"  "<<nwires_v <<"  "<<nwires_z<<endl;

  }


  // init summary tree
  art::ServiceHandle<art::TFileService> tfs;

  ftrackTree = tfs->make<TTree>("tracks","Check tracks");
  ftrackTree->Branch("EventNum", &fEventNum, "EventNum/i");
  ftrackTree->Branch("TrackId", &fTrackId, "TrackId/i");
  ftrackTree->Branch("TrackLen",   &ftrackLen,   "TrackLen/F");
  ftrackTree->Branch("Dx",  &ftrackDx,  "Dx/F");
  ftrackTree->Branch("Dy",  &ftrackDy,  "Dy/F");
  ftrackTree->Branch("Dz",  &ftrackDz,  "Dz/F");
  ftrackTree->Branch("StartX",  &ftrackStartX,  "StartX/F");
  ftrackTree->Branch("StartY",  &ftrackStartY,  "StartY/F");
  ftrackTree->Branch("StartZ",  &ftrackStartZ,  "StartZ/F");
  ftrackTree->Branch("StartTick",  &ftrackStartTick,  "StartTick/i");
  
  ftrackTree->Branch("EndX",  &ftrackEndX,  "EndX/F");
  ftrackTree->Branch("EndY",  &ftrackEndY,  "EndY/F");
  ftrackTree->Branch("EndZ",  &ftrackEndZ,  "EndZ/F");
  ftrackTree->Branch("EndTick",  &ftrackEndTick,  "EndTick/i");

  ftrackTree->Branch("theta",&ftracktheta,  "theta/F");
  ftrackTree->Branch("phi",  &ftrackphi,  "phi/F");
  ftrackTree->Branch("norm", &ftracknorm,  "norm/F");


  frecoTree  = tfs->make<TTree>("reco","Check reconstruction");
  frecoTree->Branch("X", &frecoX, "X/F");
  frecoTree->Branch("Y", &frecoY, "Y/F");
  frecoTree->Branch("Z", &frecoZ, "Z/F");
  frecoTree->Branch("EventNum", &fEventNum, "EventNum/i");
  frecoTree->Branch("TrackId",  &fTrackId, "TrackId/i");
  frecoTree->Branch("Plane",  &frecoPlane, "Plane/i");
  frecoTree->Branch("TPC", &frecoTPC, "TPC/i");
  frecoTree->Branch("Wire",  &frecoWire, "Wire/i");
  frecoTree->Branch("Channel", &frecoChannel, "Channel/i");
  frecoTree->Branch("Time", &frecoTime, "Time/i");

}

//
void pdvdana::CheckHitsAndTracks::endJob()
{
  const string myname = "pdvdana::CheckHitsAndTracks::endJob: ";
  cout<<myname<<"tracks processed total     : "<<fTotalTracks<<endl;
  cout<<myname<<"tracks processed selected  : "<<fSelectedTracks<<endl;

}

//
// check that hits associated to a track belong to the same CRP
void pdvdana::CheckHitsAndTracks::checkTrackCharacs( float Dx, float Dy, float Dz, vector<float> &TrackCharac ){
  
  double norm = sqrt(Dx*Dx+Dy*Dy+Dz*Dz);
  if (Dx>0){
    Dx = -Dx;
    Dz = -Dz;
  }
  TrackCharac[0] = Dx;
  TrackCharac[1] = Dy;
  TrackCharac[2] = Dz;
  TrackCharac[3] = norm;
  TrackCharac[4] = 0;
  TrackCharac[5] = 0;

  if (norm<1)
    return;

  // Determination of the angle
  double theta = 180 - acos(abs(Dx)/norm)*180/3.1415;
  double phi = 0;

  if(Dy>=0){
    if (Dz>=0)
      phi = atan(abs(Dz)/abs(Dy))*180/3.1415 ;
    else{
      phi = 270+atan(abs(Dy)/abs(Dz))*180/3.1415;
    }
  }
  else{
    if (Dz<0){
      phi= atan(abs(Dz)/abs(Dy))*180/3.1415+180;
    }
    else{
      phi= 90+atan(abs(Dy)/abs(Dz))*180/3.1415;
    }
  }
  TrackCharac[4] = theta;
  TrackCharac[5] = phi;
  return;
}
void pdvdana::CheckHitsAndTracks::DrawSquare(TPolyLine *Boite, float a_1, float a_2, float b_1, float b_2){    
  Boite->SetPoint(0,a_1, b_1);
  Boite->SetPoint(1,a_1, b_2);    
  Boite->SetPoint(2,a_2, b_2);
  Boite->SetPoint(3,a_2, b_1);
  Boite->SetPoint(4,a_1, b_1);
  Boite->SetLineColor(12);
  Boite->SetLineWidth(3);
  return;
}

void pdvdana::CheckHitsAndTracks::DrawCube(TPolyLine3D *Boite, float x_min, float y_min, float z_min, float x_max, float y_max, float z_max ){
  Boite->SetPoint(0,x_min, y_min, z_min);
  Boite->SetPoint(1,x_max, y_min, z_min);
  Boite->SetPoint(2,x_max, y_max, z_min);
  Boite->SetPoint(3,x_min, y_max, z_min);
  Boite->SetPoint(4,x_min, y_min, z_min);
  Boite->SetPoint(5,x_min, y_min, z_max);
  Boite->SetPoint(6,x_max, y_min, z_max);
  Boite->SetPoint(7,x_max, y_min, z_min);
  Boite->SetPoint(8,x_max, y_max, z_min);
  Boite->SetPoint(9,x_max, y_max, z_max);
  Boite->SetPoint(10,x_min, y_max, z_max);
  Boite->SetPoint(11,x_min, y_min, z_max);
  Boite->SetPoint(12,x_min, y_max, z_max);    
  Boite->SetPoint(13,x_min, y_max, z_min);
  Boite->SetPoint(14,x_min, y_max, z_max);
  Boite->SetPoint(15,x_max, y_max, z_max);
  Boite->SetPoint(16,x_max, y_min, z_max);
  Boite->SetLineColor(12);
  Boite->SetLineWidth(3);
  return;
}

void pdvdana::CheckHitsAndTracks::SetHitStyle(TPolyMarker *l_hits){
  l_hits->SetMarkerStyle(8);
  l_hits->SetMarkerSize(0.5);    
  l_hits->SetMarkerColor(kBlue);
  return;
}
unsigned pdvdana::CheckHitsAndTracks::GetWireOffset(unsigned plane_id, unsigned tpc_id){

  if (plane_id==0){
    return fOffsetWireID_u[tpc_id];
  }
  else if (plane_id==1){
    return fOffsetWireID_v[tpc_id];
  }
  else if(plane_id==2){
    return fOffsetWireID_z[tpc_id];
  }

  return 0;

}

DEFINE_ART_MODULE(pdvdana::CheckHitsAndTracks)
