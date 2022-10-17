#ifndef THINSLICESYSTEMATIC_hh
#define THINSLICESYSTEMATIC_hh
#include "fhiclcpp/ParameterSet.h"
namespace protoana {
class ThinSliceSystematic {
 public:
  ThinSliceSystematic(const fhicl::ParameterSet & pset)
    : fType(pset.get<std::string>("Name")),
      fVal(pset.get<double>("Central")),
      fCentral(pset.get<double>("Central")),
      fUpperLimit(pset.get<double>("UpperLimit")),
      fLowerLimit(pset.get<double>("LowerLimit")),
      fThrowLimit(pset.get<double>("ThrowLimit")),
      fThrowLimitUp(pset.get<double>("ThrowLimitUp")),
      fGenThrowLimit(pset.get<double>("GenThrowLimit")),
      fGenThrowLimitUp(pset.get<double>("GenThrowLimitUp")),
      fSigma(pset.get<double>("Sigma", 1.)),
      fIsG4RWCoeff(pset.get<bool>("IsG4RWCoeff", false)),
      fIsTiedG4RWCoeff(pset.get<bool>("IsTiedG4RWCoeff", false)),
      fOptions(pset.get<fhicl::ParameterSet>("Options")) {
    fName = "par_" + fType + "_syst";
    if (fIsG4RWCoeff) {
      fG4RWCoeffBranch = fOptions.get<std::string>("Branch");
    }
    if (fIsTiedG4RWCoeff) {
      fTiedG4RWCoeffBranches
          = fOptions.get<std::vector<std::string>>("Branches");
    }
  };

  ThinSliceSystematic(){};

  ~ThinSliceSystematic(){};

  const double GetValue() const {
    return fVal;
  };

  const std::string & GetG4RWCoeffBranch() const {
    return fG4RWCoeffBranch;
  };

  const std::vector<std::string> & GetTiedG4RWBranches() const {
    return fTiedG4RWCoeffBranches;
  };

  const double GetCentral() const {
    return fCentral;
  };

  const double GetUpperLimit() const {
    return fUpperLimit;
  };

  const double GetThrowLimit() const {
    return fThrowLimit;
  };
  const double GetThrowLimitUp() const {
    return fThrowLimitUp;
  };

  const double GetGenThrowLimit() const {
    return fGenThrowLimit;
  };
  const double GetGenThrowLimitUp() const {
    return fGenThrowLimitUp;
  };

  const double GetLowerLimit() const {
    return fLowerLimit;
  };

  const std::string GetName() const {
    return fName;
  };

  void SetValue(double v) {
    fVal = v;
  };

  void SetToNSigma(int n = 1) {
    SetValue(n*fSigma);
  }

  void SetToCentral() {
    SetValue(fCentral); 
  }

  template <typename T> const T GetOption(std::string name) const {
    return fOptions.get<T>(name);
  }

  const bool GetIsG4RWCoeff() const {
    return fIsG4RWCoeff;
  }

  const bool GetIsTiedG4RWCoeff() const {
    return fIsTiedG4RWCoeff;
  }
 private:

  std::string fType;
  double fVal;
  double fCentral;
  double fUpperLimit; 
  double fLowerLimit;
  double fThrowLimit; 
  double fThrowLimitUp; 
  double fGenThrowLimit; 
  double fGenThrowLimitUp; 
  double fSigma;

  bool fIsG4RWCoeff, fIsTiedG4RWCoeff;
  std::string fG4RWCoeffBranch = "";
  std::vector<std::string> fTiedG4RWCoeffBranches;

  fhicl::ParameterSet fOptions;

  std::string fName;
};
}
#endif
