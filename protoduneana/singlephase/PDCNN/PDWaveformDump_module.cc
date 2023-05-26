////////////////////////////////////////////////////////////////////////
// Class:       PDWaveformDump
// Plugin Type: analyzer (Unknown Unknown)
// File:        PDWaveformDump_module.cc
//
// Generated at Tue May  3 16:52:11 2022 by Tingjun Yang using cetskelgen
// Based on LED analysis code written by Dante Totani
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RawData/OpDetWaveform.h"

#include "TTree.h"

#include <iostream>
#include <vector>

using namespace std;

namespace pdsp {
  class PDWaveformDump;
}


class pdsp::PDWaveformDump : public art::EDAnalyzer {
public:
  explicit PDWaveformDump(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PDWaveformDump(PDWaveformDump const&) = delete;
  PDWaveformDump(PDWaveformDump&&) = delete;
  PDWaveformDump& operator=(PDWaveformDump const&) = delete;
  PDWaveformDump& operator=(PDWaveformDump&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;

private:

  // Declare member data here.
  TTree *pdtree;
  int run; // run number
  int event; // event number
  int daqch; // channel number
  vector<short> onda; // waveform
};


pdsp::PDWaveformDump::PDWaveformDump(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void pdsp::PDWaveformDump::analyze(art::Event const& e)
{
  // Get OpDetWaveform
  auto wfListHandle = e.getHandle< std::vector<raw::OpDetWaveform> >("ssprawdecoder:external");
  if (!wfListHandle){
    cout<<"wfListHandle invalid"<<endl;
    return;
  }

  run = e.run();
  event = e.id().event();
  for (auto const & wf : *wfListHandle){
    daqch = wf.ChannelNumber();
    onda.clear();
    for (unsigned short i = 0; i<wf.Waveform().size(); ++i){
      onda.push_back(wf.Waveform()[i]);
    }
    pdtree->Fill();
  }
}

void pdsp::PDWaveformDump::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  pdtree = tfs->make<TTree>("pdtree", "PD waveform");
  pdtree->Branch("run", &run);
  pdtree->Branch("event", &event);
  pdtree->Branch("daqch", &daqch);
  pdtree->Branch("onda", &onda);
}

DEFINE_ART_MODULE(pdsp::PDWaveformDump)
