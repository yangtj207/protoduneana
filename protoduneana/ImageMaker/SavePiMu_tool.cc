#include "art/Utilities/ToolMacros.h" 
#include "ImageMaker.h" 

namespace dnn{

  void saveImage(art::Event const& e, hep_hpc::hdf5::File &hdffile){
    return;
  }
}

DEFINE_ART_FUNCTION_TOOL(dnn::saveImage, "pimu")
