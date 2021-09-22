#ifndef PROTODUNE_EMPTY_EVENT_FINDER_H
#define PROTODUNE_EMPTY_EVENT_FINDER_H

#include <string>
#include <vector>

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"

namespace protoana
{
  class ProtoDUNEEmptyEventFinder
  {
  
    public:
  
    ProtoDUNEEmptyEventFinder(const fhicl::ParameterSet &pset);
  
    bool IsEmptyEvent(const art::Event &evt) const;
  
    private:
  
    // Input label for the space points. By default use the space point solver
    // as it is pseudo independent of Pandora
    std::string fSpacePointLabel;
  
    // Number of hits required in the beam region for an event to not 
    // be classed as empty
    unsigned int fNHitsThreshold;
  
    // Store the x, y and z values to define the beam region cuboid
    std::vector<float> fMinCoordsData;
    std::vector<float> fMaxCoordsData;
    std::vector<float> fMinCoordsSim;
    std::vector<float> fMaxCoordsSim;
  };
}
#endif
