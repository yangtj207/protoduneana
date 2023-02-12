#include "art/Utilities/ToolMacros.h" 
#include "ImageMaker.h" 
#include "hep_hpc/hdf5/Ntuple.hpp"
#include "TTimeStamp.h"

namespace dnn{

  inline std::array<unsigned int, 3>
  get_eid(art::Event const& e)
  {
    return { e.run(), e.subRun(), e.id().event() };
  }
  
  void saveImage(art::Event const& e, hep_hpc::hdf5::File &hdffile){
    static hep_hpc::hdf5::Ntuple<unsigned int, double> evtids(hdffile, "evtids", {{"eid",3}, "evttime"});
    
    auto event_id = get_eid(e);
    double evttime;

    art::Timestamp ts = e.time();
    if (ts.timeHigh() == 0){
      TTimeStamp tts(ts.timeLow());
      evttime = tts.AsDouble();
    }
    else{
      TTimeStamp tts(ts.timeHigh(), ts.timeLow());
      evttime = tts.AsDouble();
    }
    evtids.insert(event_id.data(), evttime);

    return;
  }
}

DEFINE_ART_FUNCTION_TOOL(dnn::saveImage, "pimu")
