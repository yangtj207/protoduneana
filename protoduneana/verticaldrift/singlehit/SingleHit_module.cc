////////////////////////////////////////////////////////////////////////
// Class:       SingleHit
// Plugin Type: analyzer (Unknown Unknown)
// File:        SingleHit_module.cc
//
// Generated at Thu Feb  1 04:19:28 2024 by Emile Lavaut using cetskelgen
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

namespace pdvdana {
  class SingleHit;
}


class pdvdana::SingleHit : public art::EDAnalyzer {
public:
  explicit SingleHit(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SingleHit(SingleHit const&) = delete;
  SingleHit(SingleHit&&) = delete;
  SingleHit& operator=(SingleHit const&) = delete;
  SingleHit& operator=(SingleHit&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
   Tree *fAnaTree;

  // Tree variables
  unsigned int fEventID;
  int fPlane;
  int fChannel;
  int fHitWidth;
  int fCutGeometry;
  int fCutHit;
  int fCoincidence;

  float fEnergy;
  float fPeakTime;
  float fSigmaPeakTime;
  float fRMS;
  float fAmplitude;
  float fSigmaAmplitude;
  float fGoodnessOfFit;
  float fIntegral;
  float fSigmaIntegral;

  float fX;
  float fY;
  float fZ;

  float fEnergyInd1;
  float fPeakTimeInd1;
  float fEnergyInd2;
  float fPeakTimeInd2;

  //Input variables
  std::string fHitLabel;
  std::string fSpacePointLabel;
  std::string fClusterLabel;
  std::string fTrackLabel;

  int fChannelWd;

  float fPeakTimeWd;
  float fCoincidenceWd;

  float fGeometryZWd;
  float fGeometryYWd;

};


pdvdana::SingleHit::SingleHit(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fHitLabel(p.get<std::string>("HitLabel")),
    fSpacePointLabel(p.get<std::string>("SpacePointLabel")),
    fClusterLabel(p.get<std::string>("ClusterLabel")),
    fTrackLabel(p.get<std::string>("TrackLabel")),

    fChannelWd(p.get<int>("ChannelWindow")),
    fPeakTimeWd(p.get<float>("PeakTimeWindow")),
    fCoincidenceWd(p.get<float>("CoincidenceWindow")),
    fGeometryZWd(p.get<float>("GeometryZWindow")),
    fGeometryYWd(p.get<float>("GeometryYWindow"))  // ,
    // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void pdvdana::SingleHit::analyze(art::Event const& e)
{
  //Set event ID
  fEventID = e.id().event();

  //Filling tree
  fAnaTree->Fill();

  // Implementation of required member function here.
}

void pdvdana::SingleHit::beginJob()
{
  // Implementation of optional member function here.
  art::ServicesHandle<art::TFileService> tfs;
  fAnaTree = tfs->make<TTree>("anaTree", "Output Tree");

  //add branches
  fAnaTree->Branch("eventID"    , &fEventID    );
  fAnaTree->Branch("plane"      , &fPlane      );
  fAnaTree->Branch("channel"    , &fChannel    );
  fAnaTree->Branch("cutGeometry", &fCutGeometry);
  fAnaTree->Branch("cutHit"     , &fCutHit     );
  fAnaTree->Branch("coincidence", &fCoincidence);

  fAnaTree->Branch("energy"           , &fEnergy        );
  fAnaTree->Branch("peakTime"         , &fPeaktime      );
  fAnaTree->Branch("hitWidth"         , &fHitWidth      );
  fAnaTree->Branch("sigmaPeakTime"    , &fSigmaPeakTime );
  fAnaTree->Branch("amplitudePeaktime", &fAmplitude     );
  fAnaTree->Branch("sigmaAmplitude"   , &fSigmaAmplitude);
  fAnaTree->Branch("rms"              , &fRMS           );
  fAnaTree->Branch("goodnessOfFit"    , &GoodnessOfFit  );
  fAnaTree->Branch("integral"         , &Integral       );
  fAnaTree->Branch("sigmaIntegral"    , &SigmaIntegral  );
  fAnaTree->Branch("x"                , &fX             );
  fAnaTree->Branch("y"                , &fY             );
  fAnaTree->Branch("z"                , &fZ             );

  fAnaTree->Branch("energyInd1"       , &fEnergyInd1    );
  fAnaTree->Branch("peakTimeInd1"     , &fPeakTimeInd1  );
  fAnaTree->Branch("energyInd2"       , &fEnergyInd2    );
  fAnaTree->Branch("peakTimeInd2"     , &fPeakTimeInd2  );

}

void pdvdana::SingleHit::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(pdvdana::SingleHit)
