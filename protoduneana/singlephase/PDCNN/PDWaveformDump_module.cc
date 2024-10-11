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
#include "canvas/Persistency/Common/FindOneP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpWaveform.h"
#include "lardataobj/RecoBase/Hit.h"

#include "TTree.h"
#include "TH1D.h"

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
  vector<short> waveform; // waveform
  vector<float> wiener;   // wiener filtered
  vector<float> impulse;  // impulse after removing responses
  TTree *hittree;
  int channel;
  float peaktime;
  float rms;
  float peak;
  float integral;
  int multiplicity;
  int localindex;
  int starttick;
  int endtick;
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
  //cout<<"run "<<run<<" event "<<event<<endl;
  auto op1ListHandle = e.getHandle< std::vector<recob::OpWaveform> >("calpd");
  auto op2ListHandle = e.getHandle< std::vector<recob::OpWaveform> >("calpd:wiener");
  //cout<<op1ListHandle->size()<<" "<<op2ListHandle->size()<<endl;
  art::FindOneP<recob::OpWaveform> assn(wfListHandle, e, "calpd");
  //cout<<assn.size()<<endl;
  art::FindOneP<recob::OpWaveform> assn2(wfListHandle, e, "calpd:wiener");
  //cout<<assn2.size()<<endl;

  for (size_t i = 0; i<(*wfListHandle).size(); ++i){
    auto &wf = (*wfListHandle)[i];
//    auto &op1 = (*op1ListHandle)[i];
//    auto &op2 = (*op2ListHandle)[i];
    //for (auto const & wf : *wfListHandle){
    daqch = wf.ChannelNumber();
    waveform.clear();
    wiener.clear();
    impulse.clear();
    //cout<<i<<" "<<assn.at(i).isAvailable()<<" "<<assn2.at(i).isAvailable()<<endl;
//    std::vector<short> temp;
//    for (unsigned short j = 0; j<wf.Waveform().size(); ++j){
//      temp.push_back(wf.Waveform()[j]);
//    }
//    sort(temp.begin(), temp.end());
//    int ped = temp[temp.size()/2];
    double baseline = 0;
    TH1D *basehelp= new TH1D("basehelp","basehelp",2000, 1300,1800);
    for(size_t j=0; j<wf.Waveform().size(); j++){
      if(j<1000){
        basehelp->Fill(wf.Waveform()[j]);
      }
    }
    int basebinmax = basehelp->GetMaximumBin();
    baseline = basehelp->GetXaxis()->GetBinCenter(basebinmax);
    basehelp->Delete();  

    for (unsigned short j = 0; j<wf.Waveform().size(); ++j){
      waveform.push_back(wf.Waveform()[j] - baseline);
//      wiener.push_back(op2.Signal()[j]);
//      impulse.push_back(op1.Signal()[j]);
      //cout<<wf.Waveform()[j]<<" "<<op1.Signal()[j]<<" "<<op2.Signal()[j]<<" "<<wf.ChannelNumber()<<" "<<op1.Channel()<<" "<<op2.Channel()<<endl;
    }
    if (assn.at(i).isAvailable()){
      auto &op = assn.at(i);
      //cout<<op.key()<<" "<<op.id().value()<<endl;
      for (unsigned short j = 0; j<op->Signal().size(); ++j){
        impulse.push_back(op->Signal()[j]);
      }
    }
    if (assn2.at(i).isAvailable()){
      auto &op = assn2.at(i);
      //cout<<op.key()<<" "<<op.id().value()<<endl;
      for (unsigned short j = 0; j<op->Signal().size(); ++j){
        wiener.push_back(op->Signal()[j]);
      }
    }
//    ip.clear();
//    if (assn.isValid()){
//      auto const &wf1 = assn.at(i);
//      cout<<"wf1="<<wf1<<endl;
//      if (!wf1) continue;
//      for (size_t j = 0; j< wf1.ref().Signal().size(); ++j){
//        ip.push_back(wf1.ref().Signal()[j]);
//      }
//    }
//    wn.clear();
//    if (assn2.isValid()){
//      auto const &wf2 = assn2.at(i);
//      if (!wf2) continue;
//      cout<<"wf2="<<wf2<<endl;
//      for (size_t j = 0; j< wf2.ref().Signal().size(); ++j){
//        wn.push_back(wf2.ref().Signal()[j]);
//      }
//    }
    pdtree->Fill();
  }

  auto hitListHandle = e.getHandle< std::vector<recob::Hit> >("gausophit");
  for (size_t i = 0; i<(*hitListHandle).size(); ++i){
    auto &hit = (*hitListHandle)[i];
    channel = hit.Channel();
    peaktime = hit.PeakTime();
    rms = hit.RMS();
    peak = hit.PeakAmplitude();
    integral = hit.Integral();
    multiplicity = hit.Multiplicity();
    localindex = hit.LocalIndex();
    starttick = hit.StartTick();
    endtick = hit.EndTick();
    hittree->Fill();
  }
}

void pdsp::PDWaveformDump::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  pdtree = tfs->make<TTree>("pdtree", "PD waveform");
  pdtree->Branch("run", &run);
  pdtree->Branch("event", &event);
  pdtree->Branch("daqch", &daqch);
  pdtree->Branch("waveform", &waveform);
  pdtree->Branch("wiener", &wiener);
  pdtree->Branch("impulse", &impulse);
  hittree = tfs->make<TTree>("hittree", "hittree");
  hittree->Branch("channel", &channel);
  hittree->Branch("peaktime", &peaktime);
  hittree->Branch("rms", &rms);
  hittree->Branch("peak", &peak);
  hittree->Branch("integral", &integral);
  hittree->Branch("multiplicity", &multiplicity);
  hittree->Branch("localindex", &localindex);
  hittree->Branch("starttick", &starttick);
  hittree->Branch("endtick", &endtick);
                 
}

DEFINE_ART_MODULE(pdsp::PDWaveformDump)
