
#ifndef __CutsetTrainer__HH
#define __CutsetTrainer__HH

#include "TFile.h"
#include "TH1.h"
#include "TGraph.h"
#include "TObject.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TTree.h"
#include "TLegend.h"
#include "TCut.h"
#include "TEventList.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStopwatch.h"
#include "Math/ProbFuncMathCore.h"

class CutsetTrainer : public TObject {
public:
  double stepSize_; //fraction of signal or background to cut in each step
  TTree* sigTree_; //tree containing signal ntuple
  TTree* bkgTree_; //tree containing background ntuple
  double endLoss_; // final loss allowed before exiting
  double tolerance_; //tolerance value for step size
  std::vector<double> maxVal_; //max values for cuts
  std::vector<double> minVal_; //min values for cuts
  std::vector<int>    cutMin_; //whether or not to cut on min
  std::vector<int>    cutMax_; //whether or not to cut on max
  std::vector<TString> variables_; //list of variables
  TString cuts_; //cut string
  TString sigID_; //cut to ID the signal
  TString bkgID_; //cut to ID the background
  int verbose_; //adds extra print outs
  TGraph* roc_; //graph of signal efficiency vs background rejection
  TGraph* func_; //graph of function values vs signal efficiency
  std::vector<double> sigEff_; //list of signal efficiencies
  std::vector<double> bkgEff_; //list of background rejections
  std::vector<TString> cutPath_; //list of cut strings used in ROC
  std::vector<double> functionVals_; //list of optimizing function values
  TString sigWeight_;
  TString bkgWeight_;
  TString function_; //function to optimize, limit, significance, or s/b
  double lum_; //luminosity scale to yield numbers
  TStopwatch* timer_;


public:

  CutsetTrainer() :
    stepSize_(0.02),
    endLoss_(0.2),
    tolerance_(0.001),
    cuts_(""),
    sigID_(""),
    bkgID_(""),
    verbose_(0),
    sigWeight_(""),
    bkgWeight_(""),
    function_("s/b"),
    lum_(1.)
  {
    timer_ = new TStopwatch();
  }

  CutsetTrainer(double stepsize) :
    CutsetTrainer()
  {
    stepSize_ = stepsize;
  }
  
  //Setters for the fields
  void SetStepSize(double  stepsize)      {stepSize_=stepsize;}
  void SetBkgTree (TTree*  bkgTree )      {bkgTree_ =bkgTree ;}
  void SetSigTree (TTree*  sigTree )      {sigTree_ =sigTree ;}
  void SetEndLoss (double  endloss )      {endLoss_ =endloss ;}
  void SetSigID   (TString sigID   )      {sigID_   =sigID   ;}
  void SetBkgID   (TString bkgID   )      {bkgID_   =bkgID   ;}
  void SetFunction(TString function)      {function_=function;}
  void SetVerbose (int     verbose )      {verbose_ =verbose ;}

  //Getters for the fields
  double  GetStepSize()      {return stepSize_;}
  TTree*  GetBkgTree ()      {return bkgTree_ ;}
  TTree*  GetSigTree ()      {return sigTree_ ;}
  double  GetEndLoss ()      {return endLoss_ ;}
  TString GetSigID   ()      {return sigID_   ;}
  TString GetBkgID   ()      {return bkgID_   ;}
  TString GetFunction()      {return function_;}
  int     GetVerbose ()      {return verbose_ ;}

  virtual int Train();
  virtual double CutStep();
  virtual double GetNumPass(TString selection, TTree* tree, TString func);
  virtual double EvaluateFunction(double nsignal, double nbackground);
  virtual int AddVariable(TString name);
  virtual int GetMaxMin();
  virtual void BuildCutString();
  virtual TCanvas* PlotROC();
  virtual TCanvas* PlotFunction();

};

#endif
