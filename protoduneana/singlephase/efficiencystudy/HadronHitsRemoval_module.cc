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
private:
  // Declare member data here.
};


pdsp::HadronHitsRemoval::HadronHitsRemoval(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void pdsp::HadronHitsRemoval::produce(art::Event& e)
{
  // Implementation of required member function here.
  // Add code to select beam tracks using Pandora information
}

DEFINE_ART_MODULE(pdsp::HadronHitsRemoval)
