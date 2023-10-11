#ifndef THINSLICEDRIVER_hh
#define THINSLICEDRIVER_hh

#include <map>
#include <vector>

#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"

#include "ThinSliceSample.h"
#include "ThinSliceDataSet.h"
#include "ThinSliceSystematic.h"
#include "ThinSliceEvent.h"

#include "fhiclcpp/ParameterSet.h"

namespace protoana {
class ThinSliceDriver {
 public:
  ThinSliceDriver(const fhicl::ParameterSet & extra_options);
  virtual ~ThinSliceDriver();

  virtual void FillMCEvents(
    TTree * tree, std::vector<ThinSliceEvent> & events,
    std::vector<ThinSliceEvent> & fake_data_events,
    int & split_val, const bool & do_split, const bool & shuffle,
    int max_entries, const bool & do_fake_data) = 0;

  virtual void BuildDataHists(
    TTree * tree, ThinSliceDataSet & data_set, double & flux,
    const std::vector<double> & beam_energy_bins,
    std::vector<double> & beam_fluxes,
    int split_val = 0) = 0;

  virtual void FillDataHistsFromSamples(
    ThinSliceDataSet & data_set,
    const std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    double & flux, std::vector<double> & fluxes_by_beam, bool fluctuate
  ) = 0;

  virtual void BuildFakeData(
    TTree * tree, const std::vector<ThinSliceEvent> & events,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    const std::map<int, bool> & signal_sample_checks,
    ThinSliceDataSet & data_set, double & flux,
    std::map<int, std::vector<double>> & sample_scales,
    std::vector<double> & beam_energy_bins,
    std::vector<double> & beam_fluxes,
    int split_val = 0, bool scale_to_data_beam_p = false) = 0;

  virtual void BuildMCSamples(
      //TTree * tree,
      const std::vector<ThinSliceEvent> & events,
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      const std::map<int, bool> & signal_sample_checks,
      std::map<int, double> & nominal_fluxes,
      std::map<int, std::vector<std::vector<double>>> & fluxes_by_sample,
      std::vector<double> & beam_energy_bins, bool use_beam_inst_P) = 0;

  virtual void RefillMCSamples(
      const std::vector<ThinSliceEvent> & events,
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      const std::map<int, bool> & signal_sample_checks,
      std::vector<double> & beam_energy_bins,
      const std::map<int, std::vector<double>> & signal_pars,
      const std::map<int, double> & flux_pars,
      const std::map<std::string, ThinSliceSystematic> & syst_pars,
      const std::map<std::string, ThinSliceSystematic> & g4rw_pars,
      bool fit_under_over, bool tie_under_over, bool use_beam_inst_P,
      bool fill_incident = false, std::map<int, TH1*> * fix_factors = 0x0) = 0;

  /*
  virtual void BuildSystSamples(
      TTree * tree,
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      const std::map<int, bool> & signal_sample_checks,
      std::vector<double> & beam_energy_bins) = 0;*/

  virtual std::pair<double, size_t> CalculateChi2(
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      ThinSliceDataSet & data_set) = 0;

  virtual void CompareSelections(
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      ThinSliceDataSet & data_set,
      TFile & output_file,
      std::vector<std::pair<int, int>> plot_style,
      bool plot_rebinned,
      bool post_fit, int nPars,
      TDirectory * plot_dir) = 0;

  virtual void GetCurrentHists(
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      ThinSliceDataSet & data_set,
      std::map<int, std::vector<TH1*>> & hists,
      bool plot_rebinned) = 0;

  virtual void GetCurrentTruthHists(
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      std::map<int, std::vector<TH1*>> & hists,
      std::map<int, std::vector<TH1*>> & inc_hists,
      std::map<int, std::vector<TH1*>> & xsec_hists,
      const std::vector<int> & incident_samples,
      const std::map<int, std::vector<double>> & signal_bins) = 0;

  virtual void PlotThrows(
    ThinSliceDataSet & data_set, std::map<int, std::vector<TH1*>> & throw_hists,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    size_t nThrows,
    std::map<int, std::vector<TH1*>> & truth_throw_hists,
    std::map<int, std::vector<TH1*>> & truth_inc_hists,
    std::map<int, std::vector<TH1*>> & truth_xsec_hists,
    std::map<int, TH1*> & best_fit_incs,
    std::map<int, TH1*> & best_fit_xsecs,
    std::map<int, TH1*> & nominal_incs,
    std::map<int, TH1*> & nominal_xsecs,
    TFile & output_file, bool plot_rebinned,
    std::map<int, std::vector<double>> * sample_scales = 0x0) = 0;

  /*virtual void PostFitThrows(
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      ThinSliceDataSet & data_set,
      TH1D & pars,
      TH2D & cov,
      TFile & output_file,
      std::vector<std::pair<int, int>> plot_style,
      bool plot_rebinned) = 0;*/

  void CompareDataMC(
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      ThinSliceDataSet & data_set,
      TFile & output_file,
      std::vector<std::pair<int, int>> plot_style,
      int nPars,
      TDirectory * plot_dir,
      bool plot_rebinned = false,
      bool post_fit = false);

  std::pair<int, int> GetColorAndStyle(
      size_t i, const std::vector<std::pair<int, int>> & plot_style);

  virtual void SetupSysts(
      const std::vector<ThinSliceEvent> & events,
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      const std::map<int, bool> & signal_sample_checks,
      std::vector<double> & beam_energy_bins,
      const std::map<std::string, ThinSliceSystematic> & pars,
      const std::map<std::string, ThinSliceSystematic> & g4rw_pars,
      TFile & output_file) = 0;
  //virtual void WrapUpSysts(TFile & output_file) = 0;

  virtual void ConstructCovariances(
    const std::vector<ThinSliceEvent> & events,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & nominal_samples,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & covariance_samples,
    const std::map<int, bool> & signal_sample_checks,
    std::vector<double> & beam_energy_bins,
    std::map<int, double> & nominal_fluxes,
    std::map<int, std::vector<std::vector<double>>> & fluxes_by_sample,
    const std::map<int, std::vector<double>> & signal_pars,
    const std::map<int, double> & flux_pars,
    const std::map<std::string, ThinSliceSystematic> & syst_pars,
    const std::map<std::string, ThinSliceSystematic> & g4rw_pars,
    bool fit_under_over, bool tie_under_over, bool use_beam_inst_P
    /*   const std::vector<ThinSliceEvent> & events,
       std::map<int, std::vector<std::vector<ThinSliceSample>>> & nominal_samples,
       std::map<int, std::vector<std::vector<ThinSliceSample>>> & covariance_samples,
       const std::map<int, bool> & signal_sample_checks,
       std::map<int, double> & nominal_fluxes,
       std::map<int, std::vector<std::vector<double>>> & fluxes_by_sample,
       std::vector<double> & beam_energy_bins, bool use_beam_inst_P*/) = 0;
 virtual void SetupExtraHists(
    ThinSliceDataSet & data_set, 
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    std::map<int, std::vector<std::vector<ThinSliceSample>>> & fake_samples) = 0;

 void SaveExtraHists(ThinSliceDataSet & data_set,
                     TDirectory * plot_dir) {
   plot_dir->cd();
   for (auto & hist : data_set.GetExtraHists()) {
     std::cout << "Writing " << hist.first << " " << hist.second << " " << hist.second->GetName() << std::endl;

     hist.second->Write();
   }
 };

 void SetStatVar(bool set) {fStatVar = set;};
 void SetFillFakeInMain(bool set) {fFillFakeDataInMain = set;};
 protected:
  fhicl::ParameterSet fExtraOptions;
  std::string fFakeDataRoutine;

  void ResetSamples(
      std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples);
  bool fStatVar = false;
  bool fFillFakeDataInMain = false;
 private:
};
}
#endif
