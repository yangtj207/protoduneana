////////////////////////////////////////////////////////////////////////
// Class:       SpatialCoincidence
// Plugin Type: analyzer (Unknown Unknown)
// File:        SpatialCoincidence_module.cc
//
// Generated at Mon Feb 19 08:30:32 2024 by Emile Lavaut using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//LArSoft
#include "larcore/CoreUtils/ServiceUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/ServicePack.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Hit.h"

// ROOT
#include "Math/ProbFunc.h"
#include "TStyle.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TPaveStats.h"

namespace pdvdana {
  class SpatialCoincidence;
}


class pdvdana::SpatialCoincidence : public art::EDAnalyzer {
public:
  explicit SpatialCoincidence(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SpatialCoincidence(SpatialCoincidence const&) = delete;
  SpatialCoincidence(SpatialCoincidence&&) = delete;
  SpatialCoincidence& operator=(SpatialCoincidence const&) = delete;
  SpatialCoincidence& operator=(SpatialCoincidence&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  bool test;
  geo::WireID wa;
  unsigned int ca   = 1;
  unsigned int tpca = 1;
  unsigned int pa   = 0;
  unsigned int w1   = 1;

  geo::WireID wb;
  unsigned int cb   = 1;
  unsigned int tpcb = 1;
  unsigned int pb   = 0;
  unsigned int w2   = 1;

  geo::Point_t point = geo::Point_t(-999,-999,-999);

  double y;
  double z;
  // Declare member data here.

  const geo::Geometry* fGeom;
  
};


pdvdana::SpatialCoincidence::SpatialCoincidence(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    test(p.get<bool>("TestFlag")),
    ca(p.get<unsigned int>("Cryo")),
    tpca(p.get<unsigned int>("TPCA")),
    pa(p.get<unsigned int>("PlaneA")),
    w1(p.get<unsigned int>("WireA")),
    cb(p.get<unsigned int>("Cryo")),  
    tpcb(p.get<unsigned int>("TPCB")),  
    pb(p.get<unsigned int>("PlaneB")),  
    w2(p.get<unsigned int>("WireB"))  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  fGeom    = &*art::ServiceHandle<geo::Geometry>();
}


void pdvdana::SpatialCoincidence::analyze(art::Event const& e)
{
  // Implementation of required member function here.

  wa = geo::WireID( ca , tpca , pa , w1 );
  wb = geo::WireID( cb , tpcb , pb , w2 ); 
  bool drap = fGeom->WireIDsIntersect( wa , wb , point);

  if ( drap ) 
  {
    std::cout<<"================================================================="<< std::endl;
    std::cout<<"there is an intersection : " << drap << std::endl;
    std::cout<<"y = "<< point.Y() <<" z = "<< point.Z() <<" x = "<< point.X() << std::endl;
    std::cout<<"================================================================="<< std::endl;
  }
}

void pdvdana::SpatialCoincidence::beginJob()
{
  // Implementation of optional member function here.
}

void pdvdana::SpatialCoincidence::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(pdvdana::SpatialCoincidence)
