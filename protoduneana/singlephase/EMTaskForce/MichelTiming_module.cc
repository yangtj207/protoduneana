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
#include "canvas/Utilities/InputTag.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

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

  art::InputTag fFlashModuleLabel;

};


pdsp::MichelTiming::MichelTiming(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fFlashModuleLabel(p.get< art::InputTag >("FlashModuleLabel"))
{
}

void pdsp::MichelTiming::analyze(art::Event const& e)
{
  auto flashListHandle = e.getHandle < std::vector < recob::OpFlash > >(fFlashModuleLabel);
  if (!flashListHandle) {
    return;
  }

  std::size_t i = 0;
  for (const recob::OpFlash & flash : *flashListHandle){
    std::cout<< i++ <<" "<<flash.Time()<<" "<<flash.TotalPE()<<std::endl;
  }
    

}

void pdsp::MichelTiming::beginJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(pdsp::MichelTiming)
