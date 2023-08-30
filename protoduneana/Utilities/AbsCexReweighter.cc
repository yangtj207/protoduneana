#include "AbsCexReweighter.hh"
#include "geant4reweight/ReweightBase/G4ReweightTraj.hh"
#include "Geant4/G4PionPlus.hh"

AbsCexReweighter::AbsCexReweighter(
    TFile * FSInput,
    const std::map<std::string, TH1D*> &FSScales,
    const fhicl::ParameterSet & material_pars,
    G4ReweightManager * rw_manager,
    TH1D * inputElasticBiasHist, bool fix)
  : G4Reweighter(FSInput, FSScales, material_pars, rw_manager,
                 {"cex", "abs", "other"},
                 inputElasticBiasHist, fix) {
  part_def = piplus->Definition();
  fInelastic = "pi+Inelastic";
  std::cout << "Part def: " << part_def << std::endl;
  SetupProcesses();
}

std::string AbsCexReweighter::GetInteractionSubtype(
    const G4ReweightTraj & theTraj) {

  //int nPi0 = 0;
  //int nPiPlus = 0;
  //int nPiMinus = 0;
  const auto & children = theTraj.GetChildren();
  for (auto * child : children) {
    int pdg = child->GetPDG();
    if (pdg == 111) {
      //++nPi0;
      continue;
    }
    else if (abs(pdg) == 211) {

      double momentum = DBL_MAX;
      if (child->GetNSteps() > 0) {
        auto * step = child->GetStep(0);
        momentum = sqrt(step->GetPreStepPx()*step->GetPreStepPx() +
                        step->GetPreStepPy()*step->GetPreStepPy() +
                        step->GetPreStepPz()*step->GetPreStepPz());
      }
      std::cout << pdg << " " << momentum << std::endl;
      if (momentum > fThreshold)
        //++nPiCharged;
    }
    else continue;
  }

  int nPi0     = theTraj.HasChild(111).size();
  int nPiPlus  = theTraj.HasChild(211).size();
  int nPiMinus = theTraj.HasChild(-211).size();

  if( (nPi0 + nPiPlus + nPiMinus) == 0){
    return "abs";
  }
  else if( (nPiPlus + nPiMinus) == 0 && nPi0 > 0 ){
    return "cex";
  }

  return "other";
}

AbsCexReweighter::~AbsCexReweighter(){}
