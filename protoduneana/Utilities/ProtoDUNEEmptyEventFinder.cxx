#include <string>
#include <vector>

#include "protoduneana/Utilities/ProtoDUNEEmptyEventFinder.h"

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "lardataobj/RecoBase/SpacePoint.h"

protoana::ProtoDUNEEmptyEventFinder::ProtoDUNEEmptyEventFinder(const fhicl::ParameterSet &pset) :
  fSpacePointLabel(pset.get<std::string>("SpacePointLabel","reco3d")),
  fNHitsThreshold(pset.get<unsigned int>("NHitsThreshold",10)),
  fMinCoordsData(pset.get<std::vector<float>>("MinCoordValuesData")),
  fMaxCoordsData(pset.get<std::vector<float>>("MaxCoordValuesData")),
  fMinCoordsSim(pset.get<std::vector<float>>("MinCoordValuesSim")),
  fMaxCoordsSim(pset.get<std::vector<float>>("MaxCoordValuesSim"))
{ 
}

bool protoana::ProtoDUNEEmptyEventFinder::IsEmptyEvent(const art::Event &evt) const
{
  // Get the space points
  auto spacePointHandle = evt.getValidHandle<std::vector<recob::SpacePoint>>(fSpacePointLabel);
  std::vector<art::Ptr<recob::SpacePoint>> spacePoints;
  if(spacePointHandle.isValid())
    art::fill_ptr_vector(spacePoints,spacePointHandle);

  const std::vector<float> minValues{evt.isRealData() ? fMinCoordsData : fMinCoordsSim};
  const std::vector<float> maxValues{evt.isRealData() ? fMaxCoordsData : fMaxCoordsSim};

  // Now count the number of space points
  unsigned int nBeamRegionSpacePoints{0};
  for (const art::Ptr<recob::SpacePoint> &sp : spacePoints)
  {
    const double *xyz{sp->XYZ()};

    if(xyz[0] < minValues.at(0) || xyz[0] > maxValues.at(0)) continue;
    if(xyz[1] < minValues.at(1) || xyz[1] > maxValues.at(1)) continue;
    if(xyz[2] < minValues.at(2) || xyz[2] > maxValues.at(2)) continue;

    ++nBeamRegionSpacePoints;

    // If we reach the threshold then there is no point continuing
    if (nBeamRegionSpacePoints >= fNHitsThreshold) return false;
  }

  // If we get here then we didn't meet the threshold
  return true;
}

