#ifndef PDSPSYSTEMATICS_hh
#define PDSPSYSTEMATICS_hh

#include <map>
#include "ThinSliceDriver.h"
#include "TFile.h"
#include "TGraph2D.h"

namespace protoana {

class PDSPSystematics {
 public:
  PDSPSystematics(
      const std::vector<ThinSliceEvent> & events,
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      const std::map<int, bool> & signal_sample_checks,
      std::vector<double> & beam_energy_bins,
      const std::map<std::string, ThinSliceSystematic> & pars,
      TFile & output_file);
  virtual ~PDSPSystematics(){};

  void SetupSyst_G4RWCoeff(
      const std::map<std::string, ThinSliceSystematic> & pars);
  double GetSystWeight_G4RWCoeff(
      const ThinSliceEvent & event,
      const std::map<std::string, ThinSliceSystematic> & pars);

  void SetupSyst_BeamShift(
      const std::map<std::string, ThinSliceSystematic> & pars,
      TFile & output_file);
  double GetSystWeight_BeamShift(
      const ThinSliceEvent & event,
      const std::map<std::string, ThinSliceSystematic> & pars);

  void SetupSyst_EDivWeight(
      const std::map<std::string, ThinSliceSystematic> & pars);
  double GetSystWeight_EDiv(
      const ThinSliceEvent & event,
      const std::map<std::string, ThinSliceSystematic> & pars,
      int ediv_selection_ID);

  void SetupSyst_EndZNoTrackWeight(
      const std::map<std::string, ThinSliceSystematic> & pars);
  double GetSystWeight_EndZNoTrack(
      const ThinSliceEvent & event,
      int signal_index,
      const std::map<std::string, ThinSliceSystematic> & pars,
      int upstream_ID, int no_track_ID);

  void SetupSyst_BeamMatch(
      const std::map<std::string, ThinSliceSystematic> & pars);
  double GetSystWeight_BeamMatch(
      const ThinSliceEvent & event,
      const std::map<std::string, ThinSliceSystematic> & pars,
      int upstream_ID, int no_track_ID);

  double GetSystWeight_UpstreamInt(
      const ThinSliceEvent & event,
      const std::map<std::string, ThinSliceSystematic> & pars,
      int upstream_ID);

  void SetupSyst_BoxBeam(
      const std::map<std::string, ThinSliceSystematic> & pars);
  double GetSystWeight_BoxBeam(
      const ThinSliceEvent & event,
      const std::map<std::string, ThinSliceSystematic> & pars,
      int beam_cut_ID);
 private:
  //G4RW Coeff
  std::map<std::string, std::string> fG4RWCoeffBranches;

  //BeamShift
  std::pair<double, double> fSystBeamShiftLimits;
  TGraph * fSystBeamShiftMeans, * fSystBeamShiftWidths;
  bool fSetupSystBeamShift = false;
  double fSystBeamShiftWeightCap;

  //EDiv
  double fEDivF, fEDivCut;

  //EndZ
  std::map<int, std::vector<double>> fEndZFractions;
  double fEndZNoTrackCut;

  //Beam Match
  std::vector<double> fBeamMatchLimits, fBeamMatchFractions;

  std::vector<std::pair<double, double>> fBoxBeamRegions;
  double fBoxBeamFraction;
};
}
#endif
