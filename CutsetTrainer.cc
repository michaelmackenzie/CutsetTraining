#include "cutset_training/CutsetTrainer.hh"

//trains the cut parameters in given step sizes to loss given
int CutsetTrainer::Train() {
  timer_->Reset();
  timer_->Start();
  GetMaxMin(); //get initial max/min values
  double cutFrac = 0.;
  int itr = 0;
  //TCanvas* c;
  while( cutFrac + tolerance_ < endLoss_) {
    cutFrac = CutStep(); // make a step
    itr++;
    if(cutFrac < 0) return -1;
    printf("Iteration %i with fractional loss %.3e\n",itr,cutFrac);
    // if(verbose_ > 0) c = PlotROC(); //plot ROC curve as it progresses
    // c->SetTitle("ROC Canvas");
    // c->Update();
  }
  Double_t cpuTime  = timer_->CpuTime();
  Double_t realTime = timer_->RealTime();
  printf("Processing time: %7.2fs CPU time %7.2fs Wall time\n",cpuTime,realTime);

  return 0;

}

double CutsetTrainer::CutStep() {

  
  double n = GetNumPass(sigID_, sigTree_,sigWeight_); //initial number

  //get max/min vectors that will change to find next step for each variable
  std::vector<double> maxVal = maxVal_;
  std::vector<double> minVal = minVal_;
  int vars = variables_.size();
  
  //get list of current cuts and add signal ID
  TString cutInitS = cuts_.Data();
  if(sigID_ != "" && cuts_ != "") cutInitS += "&&";
  cutInitS += sigID_;
  if(verbose_ > 0) printf("Initial signal cut = %s\n",cutInitS.Data());

  double nInit = GetNumPass(cutInitS, sigTree_,sigWeight_); //number after previous cuts
  
  //loop through each parameter
  for(int i = 0; i < vars; ++i) {
    //start with bounds min - max then binary search to find correct cut fraction
    double mx = minVal[i];
    double mn = maxVal[i];
    double npassed = 0;
    bool done = false;
    int tries = 0; // in case not possible, give it 20 attempts (due to discrete variables typically)
    int maxTries = 15;
    if(verbose_ > 0) printf("Currently examining var = %s\n", variables_[i].Data());
    //loop through with binary search until found or max attempts
    while(!done && tries < maxTries) {
      TString cutS = variables_[i].Data();
      cutS += "<";
      cutS += (mx+maxVal[i])/2.;
      if(cuts_ != "") cutS += "&&";
      cutS += cuts_.Data();
      if(sigID_ != "") cutS += "&&";
      cutS += sigID_.Data();
      npassed = GetNumPass(cutS, sigTree_,sigWeight_); // number of events to pass cuts
      tries++;
      //check if within tolerance of a step, else check if over or under the step and move in that direction
      if ( abs((nInit-npassed)*1./n - stepSize_) <= tolerance_) {maxVal[i] = (maxVal[i]+mx)/2.; done = true;}
      else if ((nInit-npassed)*1./n > stepSize_) mx = (maxVal[i]+mx)/2.;
      else if ((nInit-npassed)*1./n < stepSize_) maxVal[i] = (maxVal[i]+mx)/2.;
      //add a detailed printout for ones going high in attempts
      if(tries%15 == 0 && verbose_ > 2) printf("Try %i current mx %.3e  maxVal %.3e, nInit %.0f npassed %.0f done = %d, maxTries = %i\n",
                                               tries, mx, maxVal[i], nInit, npassed, done, maxTries);
      //if maxing out, take the overshoot to avoid never taking a step with discrete variables
      if (tries >= maxTries) {
        if(verbose_ > 0) printf("Maximum tries for variable %s, continuing with current cut param\n", variables_[i].Data());
        maxVal[i] = mx; done = true;} //decide whether or not to side to mx      
    }
    if(verbose_ > 1) printf("Finding max edge = %.3e took %i tries\n",mx,tries);
    done = false;
    tries = 0;
    //repeat the same process, but looking at cutting on the lower bound
    while(!done && tries < maxTries) {
      TString cutS(variables_[i].Data());
      cutS += ">";
      cutS += (mn+minVal[i])/2.;
      if(cuts_ != "") cutS += "&&";
      cutS += cuts_.Data();
      if(sigID_ != "") cutS += "&&";
      cutS += sigID_.Data();
      npassed = GetNumPass(cutS, sigTree_,sigWeight_); // number of events to pass cuts
      tries++;
      if (abs((nInit-npassed)*1./n - stepSize_) <= tolerance_) {minVal[i] = (minVal[i]+mn)/2.; done = true;}
      else if ((nInit-npassed)*1./n > stepSize_) mn = (minVal[i]+mn)/2.;
      else if ((nInit-npassed)*1./n < stepSize_) minVal[i] = (minVal[i]+mn)/2.;
      if(tries%15 == 0 && verbose_ > 2) printf("Try %i current mn %.3e  minVal %.3e, nInit %.0f npassed %.0f done = %d, maxTries = %i\n",
                                               tries, mn, minVal[i], nInit, npassed, done, maxTries);
      if (tries >= maxTries) {minVal[i] = mn; done = true;} //decide whether or not to side to mn
    }
    if(verbose_ > 1) printf("Finding min edge of %.3e took %i tries\n",mn,tries);
  }

  //next, evaluate each cut step to decide which is optimal
  double fracStep = 0.; //fractional loss in given step
  double fracBest = 0.; // nSig/nBkg is the current comparison
  int varBest = -1; //index of best variable
  int maxOrMin = 0; //tell if cut is below or above
  double nPassS   = 0; //for ROC curve
  double nPassB   = 0; //for ROC curve

  //loop through each variable
  for(int i = 0; i < vars; ++i) {

    //check the cut on the min value first
    TString cutS(variables_[i].Data());
    cutS += ">";
    cutS += minVal[i];
    if(cuts_ != "") cutS += "&&";
    cutS += cuts_; //common part of cut string

    TString cutSig(cutS); //string for signal, adding signal ID
    if(sigID_ != "") cutSig+="&&";
    cutSig += sigID_;

    double npassed = GetNumPass(cutSig, sigTree_, sigWeight_);  // number of events to pass cuts

    TString cutBkg(cutS.Data()); //string for background, adding background ID
    if(bkgID_ != "") cutBkg+="&&";
    cutBkg += bkgID_.Data();

    double npassedB = GetNumPass(cutBkg, bkgTree_,bkgWeight_);
    // double frac = (npassedB == 0) ? npassed*1.e9 : npassed / npassedB;
    double frac = EvaluateFunction(npassed, npassedB);
    if(verbose_ > 1) printf("Cut %s has function value %.4e from below\n", variables_[i].Data(), frac);

    //compare nSig/nBkg to previous ratios, save if best so far
    if(frac > fracBest) {
      fracBest = frac;
      nPassS  = npassed;
      nPassB  = npassedB;
      varBest = i;
      maxOrMin = -1;
      fracStep = (1.-npassed/n); // total fraction lost after this step
    }
    
    //next check cutting on the maximum
    cutS = variables_[i].Data();
    cutS += "<";
    cutS += maxVal[i];
    if(cuts_ != "") cutS += "&&";
    cutS += cuts_.Data();

    cutSig = TString(cutS.Data());
    if(sigID_ != "") cutSig +="&&";
    cutSig += sigID_.Data();
    npassed = GetNumPass(cutSig, sigTree_,sigWeight_);

    cutBkg = cutS.Data();
    if(bkgID_ != "") cutBkg+="&&";
    cutBkg += bkgID_.Data();
    npassedB = GetNumPass(cutBkg, bkgTree_,bkgWeight_);
    frac = EvaluateFunction(npassed, npassedB);
    if(verbose_ > 1) printf("Cut %s has function value %.4e from above\n", variables_[i].Data(), frac);

    //compare nSig/nBkg to previous ratios, save if best so far
    if(frac > fracBest) {
      fracBest = frac;
      nPassS  = npassed;
      nPassB  = npassedB;
      varBest = i;
      maxOrMin = 1;
      fracStep = 1.-npassed/n; // total fraction lost after this step
    }
  }

  //Get info for the best variable and set the bounds
  if(maxOrMin > 0) {
    maxVal_[varBest] = maxVal[varBest];
    cutMax_[varBest] = 1;
  } else if (maxOrMin < 0) {
    minVal_[varBest] = minVal[varBest];
    cutMin_[varBest] = 1;
  } else return -1.;
  BuildCutString();
  if(verbose_ > 0) printf("%s = %.4e (signal = %.3e, background = %.3e)\n", function_.Data(), fracBest, nPassS, nPassB);
  //add efficiencies to the ROC curve
  sigEff_.push_back(nPassS/n);
  double nB = GetNumPass(bkgID_, bkgTree_, bkgWeight_); //initial number
  bkgEff_.push_back(1.-nPassB/nB);
  if(verbose_ > 0) printf("Signal Efficiency = %.4e (%.3e / %.3e) Background rejection = %.4e (1 - %.3e / %.3e)\n", 
			  nPassS/n, nPassS, n, 1.-nPassB/nB, nPassB, nB );
  cutPath_.push_back(cuts_);
  functionVals_.push_back(fracBest); //record function vs cuts

  return fracStep;
}

double  CutsetTrainer::GetNumPass(TString selection, TTree* tree, TString func) {
  TFile* tmpFile = new TFile("CutsetTrainer_tmp_file.root","RECREATE");
  TString cutstring = selection;
  if(func != "") {
    if(cutstring != "") cutstring = Form("%s*(%s)",func.Data(),selection.Data());
    else cutstring = func;
  }
  TCut cut(cutstring);
  while(true) {
    TObject* o = gDirectory->Get("hist");
    if(!o) 
      break;
    else
      delete o;
  }
  tree->Draw(Form("%s>>hist", variables_[0].Data()), cut);
  TH1F* hist = (TH1F*) gDirectory->Get("hist");
  if(!hist) {
    printf("Failed to find %s weighted num pass histogram!", func.Data());
    return -1.;
  }
  hist->SetName("hist");
  double n = hist->Integral();
  n *= lum_;
  delete hist;
  delete tmpFile;
  return n;
}

double CutsetTrainer::EvaluateFunction(double nsignal, double nbackground) {
  double val = -1.;
  if(function_ == "s/b") {
    if(nbackground > 0.) val = nsignal/nbackground;
    else return 1.e9*nsignal;
  } else if(function_ == "significance") {
    if(nsignal > 0. || nbackground > 0.) val = nsignal/sqrt(nsignal + nbackground);
    else val = -1.;
  } else if(function_ == "limitApprox") {
    if(nbackground > 0.) val = nsignal/sqrt(nbackground)/1.64485;
    else val = 1.e9*nsignal;
  } else if(function_ == "limit") { //round background down, as n <= nbackground
    if(nbackground > 0.) {
      double p = 0.05; //confidence limit goal
      double tolerance = 0.001; //precision to achieve for goal confidence level
      double scale = 1.; //scale signal until achieve tolerance
      // double lastVal = -1.; //for searching
      // double lastScale = -1.;
      while(abs(val - p) > tolerance) { //guess scale factors until close to limit goal
	val = ROOT::Math::poisson_cdf((int) (nbackground), nbackground + nsignal*scale); //confidence limit at this value	
	if(verbose_ > 5) printf("Loop for limit: val = %.3e scale = %.3e\n", val, scale);
	if(abs(val-p) > tolerance) //only update if still not succeeding
	  scale *= (val/p < 4.) ? (1.-p)/(1.-val) : sqrt(nbackground)*1.64485/nsignal; //if far, start with this scale
      }
      val = 1./scale;
    }
    else val = 1.e9*nsignal;
  }  else {
    printf("Unknown optimizing function = %s!\n", function_.Data());
    val = -1.;
  }
  return val;
}

void CutsetTrainer::BuildCutString() {

  int vars = variables_.size();
  TString cuts = "";
  for(int var = 0; var < vars; ++var) {
    if(cutMin_[var]) {
      if(cuts != "") cuts += "&&";
      cuts += variables_[var];
      cuts += ">";
      cuts += minVal_[var];
    }
    if(cutMax_[var]) {
      if(cuts != "") cuts += "&&";
      cuts += variables_[var];
      cuts += "<";
      cuts += maxVal_[var];
    }
  }
  cuts_ = cuts;
  if(verbose_ > 1) printf("BuildCutString: Current cut string = %s \n",cuts_.Data());
}

int CutsetTrainer::AddVariable(TString name) {

  if(!sigTree_) return 1;
  if(!bkgTree_) return 2;
  if(!sigTree_->GetBranch(name.Data())) return 3;
  if(!bkgTree_->GetBranch(name.Data())) return 4;
  variables_.push_back(name);
  return 0;
}

int CutsetTrainer::GetMaxMin() {
  Float_t mx[variables_.size()];
  Float_t mn[variables_.size()];
  Float_t vars[variables_.size()];
  TTree* ttemp = sigTree_->CopyTree("");
  ttemp->SetName("ttemp");
  for(int i = 0; i < int(variables_.size()); ++i) {
    mx[i] = -1.e9;
    mn[i] =  1.e9;
    vars[i] = 0.;
    int status = ttemp->SetBranchAddress(variables_[i].Data(), &(vars[i]));
    if(status != 0) {
      if(verbose_ > 0) printf("GetMaxMin: Failed to find branch address %s!\n", variables_[i].Data());
      return i+1;
    }
  }
  Long64_t n = ttemp->GetEntriesFast();
  if(verbose_ > 2) printf("GetMaxMin: Signal tree has %lld entries\n", n);
  for(Long64_t j = 0; j < n; ++j) {
    ttemp->GetEntry(j);
    for(int i = 0; i < int(variables_.size()); ++i) {
      if     (vars[i] > mx[i]) mx[i] = vars[i];
      else if(vars[i] < mn[i]) mn[i] = vars[i];
      if(verbose_ > 2 && j == 0) printf("GetMaxMin: %s entry 0 = %.3e\n", variables_[i].Data(), vars[i]);
    }
  }
  for(int i = 0; i < int(variables_.size()); ++i) {
    if(verbose_ > 2) printf("GetMaxMin: %s initial range %.3e - %.3e\n", variables_[i].Data(), 
			    mn[i], mx[i]);
    maxVal_.push_back(mx[i]);
    minVal_.push_back(mn[i]);
    cutMin_.push_back(0);
    cutMax_.push_back(0);
  }
  delete ttemp;
  return 0;
}

TCanvas* CutsetTrainer::PlotROC() {
  int n = sigEff_.size();
  if(n == 0) return NULL;
  Double_t x[n],y[n]; 
  for(int i = 0; i < n; ++i) {
    x[i] = sigEff_[i];
    y[i] = bkgEff_[i];
  }
  roc_ = new TGraph(n,x,y);
  roc_->SetName("ROC_Graph");
  TCanvas* c = new TCanvas("ROC_Canvas","ROC Canvas");
  roc_->Draw();
  roc_->SetLineColor(kBlue);
  roc_->SetLineWidth(3);
  roc_->SetTitle("Rectangular Cuts ROC;Signal Efficiency;Background Rejection");
  return c;
}

TCanvas* CutsetTrainer::PlotFunction() {
  int n = sigEff_.size();
  if(n == 0) return NULL;
  Double_t x[n],y[n]; 
  for(int i = 0; i < n; ++i) {
    x[i] = sigEff_[i];
    y[i] = functionVals_[i];
  }
  func_ = new TGraph(n,x,y);
  func_->SetName("Function_Graph");
  TCanvas* c = new TCanvas("Function_Canvas","Function Canvas");
  func_->Draw();
  func_->SetLineColor(kBlue);
  func_->SetLineWidth(3);
  func_->SetTitle(Form("Optimizing %s;Signal Efficiency;%s", function_.Data(), function_.Data()));
  return c;
}
