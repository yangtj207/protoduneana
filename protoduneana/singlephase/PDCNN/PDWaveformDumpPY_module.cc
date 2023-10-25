////////////////////////////////////////////////////////////////////////
// Class:       PDWaveformDumpPY
// Plugin Type: analyzer (Unknown Unknown)
// File:        PDWaveformDumpPY_module.cc
//
// Generated at Tue Jul 12 12:39:29 2022 by Tingjun Yang using cetskelgen
// from  version .
// Save waveforms in numpy format for CNN training
//
////////////////////////////////////////////////////////////////////////
//THIS IS THE VERSION I AM EDITING
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larrecodnn/ImagePatternAlgs/Modules/c2numpy.h"
#include "lardataobj/RawData/OpDetWaveform.h"

#include "TH1F.h"
#include "TString.h"

using namespace std;

class PDWaveformDumpPY;


class PDWaveformDumpPY : public art::EDAnalyzer {
public:
  explicit PDWaveformDumpPY(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PDWaveformDumpPY(PDWaveformDumpPY const&) = delete;
  PDWaveformDumpPY(PDWaveformDumpPY&&) = delete;
  PDWaveformDumpPY& operator=(PDWaveformDumpPY const&) = delete;
  PDWaveformDumpPY& operator=(PDWaveformDumpPY&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  void TV1D_denoise(vector<float>& input, vector<float>& output, const int width, const float lambda);

  std::string fOutNumpyFileName;

  c2numpy_writer npywriter;

  int nsigwfs;
  int nnoisewfs;

};


PDWaveformDumpPY::PDWaveformDumpPY(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fOutNumpyFileName(p.get<std::string>("OutNumpyFileName"))
{
}

void PDWaveformDumpPY::analyze(art::Event const& e)
{
  // Get OpDetWaveform
  auto wfListHandle = e.getHandle< std::vector<raw::OpDetWaveform> >("ssprawdecoder:external");
  if (!wfListHandle){
    cout<<"wfListHandle invalid"<<endl;
    return;
  }

  std::vector<float> waveform(2000);    //raw waveform
  std::vector<float> waveformma(2000);  //moving average
  std::vector<float> waveformden(2000); //denoising
  std::vector<float> input(2000);
  std::vector<float> output(2000);

  unsigned int mobileAVG = 4;

  for (auto const & wf : *wfListHandle){

    fill(waveform.begin(), waveform.end(), 0);
    fill(waveformma.begin(), waveformma.end(), 0);
    fill(waveformden.begin(), waveformden.end(), 0);
    fill(input.begin(), input.end(), 0);
    fill(output.begin(), output.end(), 0);

    int daqch = wf.ChannelNumber();

    if ((daqch>=132 && daqch<=143) || (daqch>=96 && daqch<=107)){

      for (unsigned short i = 0; i<wf.Waveform().size(); ++i){
        waveform[i] = wf.Waveform()[i];
      }

      //========= Moving average  ============================================================
       float avg=0.;
       int c0=0, c1=0,  c2=0;
       for(size_t n=0; n<waveform.size() ; ++n){
       if(n>=mobileAVG && n<waveform.size()-mobileAVG){
          for(size_t k=n-mobileAVG; k<=n+mobileAVG; ++k){
            avg=avg+waveform[k];
            ++c0;
          }
          avg=avg/c0;
          waveformma[n]=avg;
          avg = 0;
          c0 = 0;
        }
        else{
          if(n<mobileAVG){
            for(size_t k=0; k<=n+mobileAVG; ++k){
              avg=avg+waveform[k];
              c1=c1+1;
            }
            avg=avg/c1;
            waveformma[n]=avg;
            avg = 0;
            c1 = 0;
          }
          else if(n>=waveform.size()-mobileAVG){
            for(size_t k=n-mobileAVG; k<waveform.size(); ++k){
              avg=avg+waveform[k];
              ++c2;
            }
            avg=avg/c2;
            waveformma[n]=avg;
            avg = 0;
            c2 = 0;
          }
        }
      }

      //========= Denoising ============================================================
    
       for(size_t j=0; j<waveform.size(); ++j){
	 input[j]=waveformma[j];
	 output[j]=input[j];
       }
       TV1D_denoise(input,output,2000,10);
       for(size_t j=0; j<waveform.size(); ++j){
         waveformden[j]=output[j];
       }
      
      //===========================BASELINE Histo========================================
      float base=0;
      TH1F *basehelp= new TH1F("basehelp","basehelp",2000, 1300,1800);
      for(size_t j=0; j<waveform.size(); j++){
        if(j<1000){
          basehelp->Fill(waveformden[j]);
        }
      }
      int basebinmax = basehelp->GetMaximumBin();
      base = basehelp->GetXaxis()->GetBinCenter(basebinmax);
      basehelp->Delete();

      // float sum = 0; // unused
      for(size_t j = 0; j < waveform.size(); ++j){
        waveform[j] = waveform[j]-base;
        waveformma[j] = waveformma[j]-base;
        waveformden[j] = waveformden[j]-base;
        // if (j>=1000 && j<1500) sum += waveformden[j]; // unused
      }
    
      c2numpy_uint32(&npywriter, e.id().run());
      c2numpy_uint32(&npywriter, e.id().subRun());
      c2numpy_uint32(&npywriter, e.id().event());
      c2numpy_uint8(&npywriter, daqch);
      for (size_t j = 0; j<waveform.size(); ++j){//changed to waveformma but don't think this should make difference
        c2numpy_float64(&npywriter, waveform[j]);//changed to waveform
      }
    }
  }    
}

void PDWaveformDumpPY::TV1D_denoise(vector<float>& input, vector<float>& output, const int width, const float lambda){
  if (width>0) {                /*to avoid invalid memory access to input[0]*/
    int k=0, k0=0;            /*k: current sample location, k0: beginning of current segment*/
    float umin=lambda, umax=-lambda;    /*u is the dual variable*/
    float vmin=input[0]-lambda, vmax=input[0]+lambda;    /*bounds for the segment's value*/
    int kplus=0, kminus=0;     /*last positions where umax=-lambda, umin=lambda, respectively*/
    const float twolambda=2.0*lambda;    /*auxiliary variable*/
    const float minlambda=-lambda;        /*auxiliary variable*/
    for (;;) {                /*simple loop, the exit test is inside*/
      //cout<<"k = "<<k<<endl;
      if (k>=width) return;
      while (k==width-1) {    /*we use the right boundary condition*/
        if (umin<0.0) {            /*vmin is too high -> negative jump necessary*/
          do{
            //cout<<"k0 = "<<k0<<" kminus = "<<kminus<<" vmin = "<<vmin<<endl;
            output[k0++]=vmin;
          }while (k0<=kminus);
          umax=(vmin=input[kminus=k=k0])+(umin=lambda)-vmax;
        } else if (umax>0.0) {    /*vmax is too low -> positive jump necessary*/
          do{
            //cout<<"k0 = "<<k0<<" kplus = "<<kplus<<" vmax = "<<vmax<<endl;            
            output[k0++]=vmax;
          }while (k0<=kplus);
          umin=(vmax=input[kplus=k=k0])+(umax=minlambda)-vmin;
        } else {
          vmin+=umin/(k-k0+1);
          do{
            //cout<<"k0 = "<<k0<<" k = "<<k<<" vmin = "<<vmin<<endl;                        
            output[k0++]=vmin;
          }while(k0<=k);
          return;
        }
      }
      if ((umin+=input[k+1]-vmin)<minlambda) {        /*negative jump necessary*/
        do{
          //cout<<"negative jump, k0 = "<<k0<<" vmin = "<<vmin<<endl;
          output[k0++]=vmin;
        }while (k0<=kminus);
        vmax=(vmin=input[kplus=kminus=k=k0])+twolambda;
        umin=lambda; umax=minlambda;
      } else if ((umax+=input[k+1]-vmax)>lambda) {    /*positive jump necessary*/
        do{
          //cout<<"positive jump, k0 = "<<k0<<" vmax = "<<vmax<<endl;
          output[k0++]=vmax;
        }while (k0<=kplus);
        vmin=(vmax=input[kplus=kminus=k=k0])-twolambda;
        umin=lambda; umax=minlambda;
      } else {     /*no jump necessary, we continue*/
        k++;
        if (umin>=lambda) {        /*update of vmin*/
          vmin+=(umin-lambda)/((kminus=k)-k0+1);
          umin=lambda;
        }
        if (umax<=minlambda) {    /*update of vmax*/
          vmax+=(umax+lambda)/((kplus=k)-k0+1);
          umax=minlambda;
        }
      }
    }
  }
}

void PDWaveformDumpPY::beginJob()
{

  c2numpy_init(&npywriter, fOutNumpyFileName, 50000);
  c2numpy_addcolumn(&npywriter, "run", C2NUMPY_UINT32);
  c2numpy_addcolumn(&npywriter, "subrun", C2NUMPY_UINT32);
  c2numpy_addcolumn(&npywriter, "evt", C2NUMPY_UINT32);
  c2numpy_addcolumn(&npywriter, "channel", C2NUMPY_UINT8);    
  for (int i = 0; i<2000; ++i){
    c2numpy_addcolumn(&npywriter, Form("x%d",i), C2NUMPY_FLOAT64);
  }
  nsigwfs = 0;
  nnoisewfs = 0;
}

void PDWaveformDumpPY::endJob()
{
  c2numpy_close(&npywriter);

  std::cout<<"Saved "<<nsigwfs<<" signal waveforms and "<<nnoisewfs<<" noise waveforms."<<std::endl;
}

DEFINE_ART_MODULE(PDWaveformDumpPY)
