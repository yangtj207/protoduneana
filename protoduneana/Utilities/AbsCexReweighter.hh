#ifndef AbsCexReweighter_h
#define AbsCexReweighter_h

#include "geant4reweight/ReweightBase/G4Reweighter.hh"

class G4PionPlus;

class AbsCexReweighter : public G4Reweighter {
  public:

    AbsCexReweighter(TFile *, const std::map<std::string, TH1D*> &,
                       const fhicl::ParameterSet &,
                       G4ReweightManager * rw_manager,
                       TH1D * inputElasticBiasHist = nullptr, bool fix = false);
    virtual ~AbsCexReweighter();
    std::string GetInteractionSubtype(const G4ReweightTraj &) override;
    void SetThreshold(double v) {fThreshold = v;};

  protected:
    G4PionPlus * piplus;
    double fThreshold = 150.; //150 MeV/cm threshold
};

#endif
