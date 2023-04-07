////////////////////////////////////////////////////////////////////////
// Class:       PlotSelectedTracks
// Plugin Type: analyzer (Unknown Unknown)
// File:        PlotSelectedTracks_module.cc
//
// Generated at Mon Jan 30 17:27:10 2023 by Yoann Kermaidic using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art_root_io/TFileService.h"

//LArSoft
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Track.h"

// ROOT
#include <TStyle.h>
#include "TH2D.h"
//#include "TGraph.h"
//#include "TMultiGraph.h"

using std::vector;
using std::string;

namespace pddpana {
  class PlotSelectedTracks;
}


class pddpana::PlotSelectedTracks : public art::EDAnalyzer {
public:
  explicit PlotSelectedTracks(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PlotSelectedTracks(PlotSelectedTracks const&) = delete;
  PlotSelectedTracks(PlotSelectedTracks&&) = delete;
  PlotSelectedTracks& operator=(PlotSelectedTracks const&) = delete;
  PlotSelectedTracks& operator=(PlotSelectedTracks&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  float    degTOrad = 3.14159/180.;

private:

  // Declare member data here.
  int      fLogLevel;
  string   fTrackModuleLabel;
  string   fTrackType;
  float    fFiducialCut;
  float    fTrackCurvature;
  float    fTrackDelta;
  float    fTrackDeltaStop;
  float    fTrackLenMin;
  float    fTrackLenMax;
  float    fTrackThetaMin;
  float    fTrackThetaMax;

  unsigned fTotalTracks;
  unsigned fEventNum;



  // cm - NP02 active anode volume
  float x_min = -300;
  float x_max = 300;
  float y_min = -100;
  float y_max = 0;
  float z_min = 300;
  float z_max = 400;
  float l_min = 0;
  float l_max = sqrt(pow(x_max-x_min,2)+pow(y_max-y_min,2)+pow(z_max-z_min,2));
 
  // cm - Reco uncertainty margin
  float eps   = 5;

  // cm - hit spacing threshold (holes in mis-reconstructed tracks)
  float hs_thrs = 50;

  // Store valid track start and end
  struct Edge{
    float x;
    float y;
  };

  // detector geometry
  const geo::Geometry* fGeom;

  // Store good and bad tracks hits into TH2D
  TH2D* h0;
  TH2D* h1;

  void     GetAngles(const recob::Track& track, float &theta, float &phi);
};


pddpana::PlotSelectedTracks::PlotSelectedTracks(fhicl::ParameterSet const& p)
  : EDAnalyzer{p} ,
  fLogLevel( p.get< int >("LogLevel") ),
  fTrackModuleLabel( p.get< std::string  >("TrackModuleLabel") ),
  fTrackType( p.get< std::string  >("TrackType") ),
  fFiducialCut(   p.get< float  >("FiducialCut") ),
  fTrackCurvature(p.get< float  >("TrackCurvature") ), // cm - hit to straight track spacing threshold
  fTrackDelta(    p.get< float  >("TrackDelta") ),
  fTrackDeltaStop(p.get< float  >("TrackDeltaStop") ),
  fTrackLenMin(   p.get< float  >("TrackLenMin") ),
  fTrackLenMax(   p.get< float  >("TrackLenMax") ),
  fTrackThetaMin( p.get< float  >("TrackThetaMin") ),
  fTrackThetaMax( p.get< float  >("TrackThetaMax") )
 
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  fGeom    = &*art::ServiceHandle<geo::Geometry>();
}

void pddpana::PlotSelectedTracks::analyze(art::Event const& e)
{
  bool point_set = false;
  Edge bot_edge;
  Edge top_edge;

  int hitout = 0;         // # of hits outside specs
  int isgood = 0;         // Flag good = 0 or bad = 1 tracks
  int ntrks[3] = {0,0,0};   // Number of tracks taken into account in the analysis

  // get tracks
  auto Tracks   = e.getValidHandle<vector<recob::Track>>(fTrackModuleLabel);
  fTotalTracks += Tracks->size();  
  
  if( fLogLevel >= 1 ){
    std::cout << "The event contains "<< Tracks->size() <<" tracks\n";
  }
  
  //
  fEventNum = e.id().event();

  // loop over tracks
  for (unsigned trk = 0; trk < Tracks->size(); ++trk) {
    const recob::Track& track = Tracks->at(trk);

    bot_edge.x = 0;
    bot_edge.y = 0;
    top_edge.x = 0;
    top_edge.y = 0;

    isgood = 0;
    hitout = 0;
    point_set = false;

    float theta = 0, phi = 0;
    GetAngles(track,theta,phi);

    if( fLogLevel >= 1 ){
      std::cout << "  #" << trk << " length: "<< track.Length(0) <<" cm / theta: " << theta << " deg \n";
    }

    // Filter out tracks that are shorter/longer than specified threshold
    if(track.Length(0) < fTrackLenMin || track.Length(0) > fTrackLenMax){ntrks[2]++; continue;}
      
    bool is_out = false;
    if(fTrackType == "anode"){
      // Filter out tracks that are starting outside fiducial volume
      for(size_t hit=0; hit<track.NPoints(); ++hit){
        if(track.LocationAtPoint(hit).X() == -999) continue;
        if(track.LocationAtPoint(hit).Y() - y_min < fFiducialCut || track.LocationAtPoint(hit).Y() - y_min > 100 - fFiducialCut) is_out = true;
        if(track.LocationAtPoint(hit).Z() - z_min < fFiducialCut || track.LocationAtPoint(hit).Z() - z_min > 100 - fFiducialCut) is_out = true;
        if(theta < fTrackThetaMin                              || theta > fTrackThetaMax)                                    is_out = true;
        break;
      }
    }
    else if(fTrackType == "cathode"){
      // Filter out tracks that are starting outside fiducial volume
      for(size_t hit=track.NPoints()-1; hit>0; --hit){
        if(track.LocationAtPoint(hit).X() == -999) continue;
        if(track.LocationAtPoint(hit).Y() - y_min < fFiducialCut || track.LocationAtPoint(hit).Y() - y_min > 100 - fFiducialCut) is_out = true;
        if(track.LocationAtPoint(hit).Z() - z_min < fFiducialCut || track.LocationAtPoint(hit).Z() - z_min > 100 - fFiducialCut) is_out = true;
        if(theta < fTrackThetaMin                              || theta > fTrackThetaMax                                   ) is_out = true;
        break;
      }
    }
    else if(fTrackType == "michel"){
      // Filter out tracks that are starting outside fiducial volume
      for(size_t hit=0; hit<track.NPoints(); ++hit){
        if(track.LocationAtPoint(hit).X() == -999) continue;
        if(track.LocationAtPoint(hit).Y() - y_min < fFiducialCut || track.LocationAtPoint(hit).Y() - y_min > 100 - fFiducialCut) is_out = true;
        if(track.LocationAtPoint(hit).Z() - z_min < fFiducialCut || track.LocationAtPoint(hit).Z() - z_min > 100 - fFiducialCut) is_out = true;
        break;
      }
      if(!is_out){
        // Filter out tracks that do not stop in the fiducial volume
        for(size_t hit=track.NPoints()-1; hit>0; --hit){
          if(track.LocationAtPoint(hit).X() == -999) continue;
          if(track.LocationAtPoint(hit).Y() - y_min < fFiducialCut || track.LocationAtPoint(hit).Y() - y_min > 100 - fFiducialCut) is_out = true;
          if(track.LocationAtPoint(hit).Z() - z_min < fFiducialCut || track.LocationAtPoint(hit).Z() - z_min > 100 - fFiducialCut) is_out = true;
          break;
        }
      }
    }

    if(is_out){ ntrks[2]++; continue;}

    if(track.Length(0) < l_min-eps || track.Length(0)  > l_max+eps){
      std::cout << "  BAD: Length [" << l_min-eps << "," << l_max+eps << "] -> " <<  track.Length(0) << std::endl;
      isgood = 1;
    }

    vector<bool> hitskip(track.NPoints(),false);

    // Sort tracks by goodness fTrackType defined by criteria below
    for(size_t hit=0; hit<track.NPoints(); ++hit){
      if(isgood)                                          break;
      if(track.LocationAtPoint(hit).X() == -999) {hitskip[hit] = true; continue;}

      if(!point_set){
        bot_edge.x = track.LocationAtPoint(hit).X() - x_min;
        bot_edge.y = sqrt(pow(track.LocationAtPoint(hit).Y() - y_min,2) + pow(track.LocationAtPoint(hit).Z() - z_min,2));
        point_set = true;
      }

      top_edge.x = track.LocationAtPoint(hit).X() - x_min;
      top_edge.y = sqrt(pow(track.LocationAtPoint(hit).Y() - y_min,2) + pow(track.LocationAtPoint(hit).Z() - z_min,2));

      float hit_spacing = 0;
      if(hit>0) hit_spacing = sqrt(pow(track.LocationAtPoint(hit).X()-track.LocationAtPoint(hit-1).X(),2)
                                 + pow(track.LocationAtPoint(hit).Y()-track.LocationAtPoint(hit-1).Y(),2)
                                 + pow(track.LocationAtPoint(hit).Z()-track.LocationAtPoint(hit-1).Z(),2));

      if(track.LocationAtPoint(hit).X() < x_min ||
         track.LocationAtPoint(hit).X() > x_max ||
         track.LocationAtPoint(hit).Y() < y_min-eps ||
         track.LocationAtPoint(hit).Y() > y_max+eps ||
         track.LocationAtPoint(hit).Z() < z_min-eps ||
         track.LocationAtPoint(hit).Z() > z_max+eps) {hitskip[hit] = true; hitout++;}

      if(hitout > 5 || hit_spacing > hs_thrs) {
        std::cout << "  BAD: " << hit << " " << hitout << " "
                          << track.LocationAtPoint(hit).X() << " "
                          << track.LocationAtPoint(hit).Y() << " "
                          << track.LocationAtPoint(hit).Z() << std::endl;
        isgood = 1;
        break;
      }
    }

    // Compute average track slope from the two track extremum points
    // later used to compute the minimal distance between hit point and line from extremums 

    if(top_edge.x-bot_edge.x == 0) isgood = 1;
    else if(!isgood){
      float a = 666.; a = (top_edge.y-bot_edge.y)/(top_edge.x-bot_edge.x);
      float b = 666.; b = bot_edge.y - a*bot_edge.x;

      for(size_t hit=0; hit<track.NPoints(); ++hit){
        if(hitskip[hit]) continue;

        float xpt = track.LocationAtPoint(hit).X() - x_min;
        float yzpt= sqrt(pow(track.LocationAtPoint(hit).Y() - y_min,2) + pow(track.LocationAtPoint(hit).Z() - z_min,2));

        float xx  = 1./(a*a+1)*(xpt+a*yzpt-a*b);
        float yy  = 1./a*(xpt-xx) + yzpt;

        if(sqrt(pow(xx-xpt,2) + pow(yy-yzpt,2)) > fTrackCurvature) {
          std::cout << "  BAD: " << hit << " distance: " << sqrt(pow(xx-xpt,2) + pow(yy-yzpt,2)) << " " << xpt << " " << yzpt << std::endl;
          isgood = 1;
          break;
        }
      }
    }

    if(fLogLevel >= 1)
      std::cout << "Event: #" << e.id().event() << ": track length: " << track.Length(0) << " cm / theta: " << track.Theta()/degTOrad << " deg found -> " << isgood << " status" << std::endl;

    // Fill shifted tracks points into TGraphs 
    for(size_t hit=0; hit<track.NPoints(); ++hit){
      if(hitskip[hit]) continue;

      float xpt = track.LocationAtPoint(hit).X() - x_min;
      float yzpt= sqrt(pow(track.LocationAtPoint(hit).Y() - y_min,2) + pow(track.LocationAtPoint(hit).Z() - z_min,2));

      if(fLogLevel > 1) std::cout << "    #hit: " << hit << " | " << yzpt << " / " << xpt << std::endl;
      if(isgood) h1->Fill(yzpt,xpt);
      else       h0->Fill(yzpt,xpt);
    }
    ntrks[isgood]++;
  }     
 
  if(fLogLevel >= 1){
    std::cout << "Total # of " << fTrackType << " tracks: " << ntrks[0]+ntrks[1] << std::endl;
    std::cout << "  - good: " << ntrks[0] << std::endl;
    std::cout << "  - bad:  " << ntrks[1] << std::endl;
    std::cout << "  - out:  " << ntrks[2] << std::endl;
    std::cout << " " << std::endl;
  }
}

void pddpana::PlotSelectedTracks::beginJob()
{
  if(fLogLevel >= 1){
    std::cout << " TrackMinLen:   " << fTrackLenMin << std::endl;
    std::cout << " TrackMaxLen:   " << fTrackLenMax << std::endl;
    std::cout << " TrackThetaMin: " << fTrackThetaMin << std::endl;
    std::cout << " TrackThetaMax: " << fTrackThetaMax << std::endl;
  }

  art::ServiceHandle<art::TFileService> tfs;
  h0 = tfs->make<TH2D>("h0","Good tracks",200,0,200,600,0,600);
  h1 = tfs->make<TH2D>("h1","Bad tracks" ,200,0,200,600,0,600);
}

void pddpana::PlotSelectedTracks::endJob()
{
  std::cout << "Done " << std::endl;
  std::cout << " " << std::endl;
}

void pddpana::PlotSelectedTracks::GetAngles(const recob::Track& track, float &theta, float &phi){

  float x_fst=-1, y_fst=-1, z_fst=-1;
  float x_lst=-1, y_lst=-1, z_lst=-1;

  for(size_t hit=0; hit<track.NPoints(); ++hit){
    if(track.LocationAtPoint(hit).X() == -999) continue;
    if(x_fst == -1){
      x_fst = track.LocationAtPoint(hit).X();
      y_fst = track.LocationAtPoint(hit).Y();
      z_fst = track.LocationAtPoint(hit).Z();
    }
    x_lst = track.LocationAtPoint(hit).X();
    y_lst = track.LocationAtPoint(hit).Y();
    z_lst = track.LocationAtPoint(hit).Z();
  }

  phi   = atan((y_lst-y_fst)/(z_lst-z_fst)) / degTOrad;
  theta = atan((x_fst-x_lst)/sqrt(pow(y_lst-y_fst,2)+pow(z_lst-z_fst,2))) / degTOrad;
}

DEFINE_ART_MODULE(pddpana::PlotSelectedTracks)
