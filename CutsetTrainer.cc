#include "cutset_training/CutsetTrainer.hh"

//trains the cut parameters in given step sizes to loss given
int CutsetTrainer::Train() {
  timer_->Reset();
  timer_->Start();
  significances_.verbose_ = verbose_;

  bkgTreeCut_ = bkgTree_; //->CloneTree();
  sigTreeCut_ = sigTree_; //->CloneTree();

  auto cachesize = 100000000U; //100MB
  sigTreeCut_->SetCacheSize(cachesize);
  sigTreeCut_->AddBranchToCache("*",true);
  sigTreeCut_->LoadBaskets();
  bkgTreeCut_->SetCacheSize(cachesize);
  bkgTreeCut_->AddBranchToCache("*",true);
  bkgTreeCut_->LoadBaskets();
  if(verbose_ > 0) printf("--- %s: Trees loaded into memory!\n", __func__);

  // UpdateTrees(); //update the event lists for the trees with the IDs

  GetMaxMin(); //get initial max/min values
  initSignal_ = GetNumPass(sigID_, sigTree_,sigWeight_);
  initBackground_ = GetNumPass(bkgID_, bkgTree_, bkgWeight_); //initial number
  if(verbose_ > 0) printf("--- %s: Initial signal = %.3e and background = %.3e\n",
                          __func__, initSignal_, initBackground_);
  if(initBackground_ < 0.) {
    printf("!!! %s: Background < 0! Exiting...\n", __func__);
    return 10;
  }

  double cutFrac = 0.;
  int itr = 0;
  while(cutFrac + tolerance_ < endLoss_) {
    cutFrac = CutStep(); // make a step
    itr++;
    if(cutFrac < 0) return -1;
    printf("--- %s: Iteration %i with fractional loss %.3e\n",__func__, itr,cutFrac);
  }
  Double_t cpuTime  = timer_->CpuTime();
  Double_t realTime = timer_->RealTime();
  printf("--- Processing time: %7.2fs CPU time %7.2fs Wall time\n",cpuTime,realTime);

  return 0;

}

double CutsetTrainer::CutStep() {


  //get max/min vectors that will change to find next step for each variable
  std::vector<double> maxVal = maxVal_;
  std::vector<double> minVal = minVal_;
  int vars = variables_.size();

  //get list of current cuts and add signal ID
  TString cutInitS = cuts_.Data();
  if(sigID_ != "" && cuts_ != "") cutInitS += "&&";
  cutInitS += sigID_;
  if(verbose_ > 0) printf("--- %s: Initial signal cut = %s\n",__func__,cutInitS.Data());

  double nInit = GetNumPass(cutInitS, sigTreeCut_,sigWeight_); //number after previous cuts

  /////////////////////////////
  // Find the next cut value //
  /////////////////////////////

  TStopwatch* timer = (verbose_ > 1) ? new TStopwatch() : 0;

  //loop through each parameter
  for(int i = 0; i < vars; ++i) {
    if(verbose_ > 1) timer->Start();

    //start with bounds min - max then binary search to find correct cut fraction
    double mx = minVal[i]; //running maximum, start at min so true max + running max / 2 = middle
    double mn = maxVal[i]; //same as maximum, just flipped
    double npassed = 0;
    bool done = false;
    int tries = 0; // in case not possible, give it 20 attempts (due to discrete variables typically)
    int maxTries = 15;
    if(verbose_ > 0) printf("--- %s: Currently examining var = %s\n", __func__, variables_[i].Data());
    //loop through with binary search until found or max attempts
    while(!done && tries < maxTries && cutSide_[i] >= 0) {
      TString cutS = variables_[i].Data();
      cutS += "<";
      cutS += (mx+maxVal[i])/2.;
      if(cuts_ != "") cutS += "&&";
      cutS += cuts_.Data();
      if(sigID_ != "") cutS += "&&";
      cutS += sigID_.Data();
      npassed = GetNumPass(cutS, sigTreeCut_,sigWeight_); // number of events to pass cuts
      tries++;
      //check if within tolerance of a step, else check if over or under the step and move in that direction
      if ( abs((nInit-npassed)*1./initSignal_ - stepSize_) <= tolerance_) {maxVal[i] = (maxVal[i]+mx)/2.; done = true;}
      else if ((nInit-npassed)*1./initSignal_ > stepSize_) mx = (maxVal[i]+mx)/2.;
      else if ((nInit-npassed)*1./initSignal_ < stepSize_) maxVal[i] = (maxVal[i]+mx)/2.;
      //add a detailed printout for ones going high in attempts
      if(tries%15 == 0 && verbose_ > 2) printf("--- %s: Try %i current mx %.3e  maxVal %.3e, nInit %.0f npassed %.0f done = %d, maxTries = %i\n",
                                               __func__, tries, mx, maxVal[i], nInit, npassed, done, maxTries);
      //if maxing out, take the overshoot to avoid never taking a step with discrete variables
      if (tries >= maxTries) {
        if(verbose_ > 0) printf("!!! %s: Maximum tries for variable %s, continuing with current cut param\n", __func__, variables_[i].Data());
        maxVal[i] = mx; done = true;
      } //decide whether or not to side to mx
    }
    if(verbose_ > 1 && cutSide_[i] >= 0) {printf("--- %s: Finding max edge = %.3e took %.1fs and %i tries\n",__func__,mx,timer->RealTime(), tries); timer->Reset(); timer->Start();}
    done = false;
    tries = 0;
    //repeat the same process, but looking at cutting on the lower bound
    while(!done && tries < maxTries && cutSide_[i] <= 0) {
      TString cutS(variables_[i].Data());
      cutS += ">";
      cutS += (mn+minVal[i])/2.;
      if(cuts_ != "") cutS += "&&";
      cutS += cuts_.Data();
      if(sigID_ != "") cutS += "&&";
      cutS += sigID_.Data();
      npassed = GetNumPass(cutS, sigTreeCut_,sigWeight_); // number of events to pass cuts
      tries++;
      if (abs((nInit-npassed)*1./initSignal_ - stepSize_) <= tolerance_) {minVal[i] = (minVal[i]+mn)/2.; done = true;}
      else if ((nInit-npassed)*1./initSignal_ > stepSize_) mn = (minVal[i]+mn)/2.;
      else if ((nInit-npassed)*1./initSignal_ < stepSize_) minVal[i] = (minVal[i]+mn)/2.;
      if(tries%15 == 0 && verbose_ > 2) printf("--- %s: Try %i current mn %.3e  minVal %.3e, nInit %.0f npassed %.0f done = %d, maxTries = %i\n",
                                               __func__, tries, mn, minVal[i], nInit, npassed, done, maxTries);
      if (tries >= maxTries) {
        if(verbose_ > 0) printf("!!! %s: Maximum tries for variable %s, continuing with current cut param\n", __func__, variables_[i].Data());
        minVal[i] = mn; done = true;
      } //decide whether or not to side to mn
    }
    if(verbose_ > 1 && cutSide_[i] <= 0) {printf("--- %s: Finding min edge of %.3e took %.1fs and %i tries\n",__func__,mn,timer->RealTime(),tries); timer->Reset();}
  }

  //next, evaluate each cut step to decide which is optimal
  double fracStep = 0.; //fractional loss in given step
  double fracBest = 0.; // function value that is the current comparison
  int varBest = -1; //index of best variable
  int maxOrMin = 0; //tell if cut is below or above
  double nPassS   = 0; //for ROC curve
  double nPassB   = 0; //for ROC curve

  //loop through each variable
  for(int i = 0; i < vars; ++i) {
    if(verbose_ > 1) timer->Start();

    //check the cut on the min value first
    TString cutS(variables_[i].Data());
    cutS += ">";
    cutS += minVal[i];
    if(cuts_ != "") cutS += "&&";
    cutS += cuts_; //common part of cut string

    TString cutSig(cutS); //string for signal, adding signal ID
    if(sigID_ != "") cutSig+="&&";
    cutSig += sigID_;

    double npassed = (cutSide_[i] <= 0) ? GetNumPass(cutSig, sigTreeCut_, sigWeight_) : -1.;  // number of events to pass cuts

    TString cutBkg(cutS.Data()); //string for background, adding background ID
    if(bkgID_ != "") cutBkg+="&&";
    cutBkg += bkgID_.Data();

    double npassedB = GetNumPass(cutBkg, bkgTreeCut_,bkgWeight_);
    double frac = (cutSide_[i] <= 0) ? EvaluateFunction(npassed, npassedB, function_) : -1.;
    if(verbose_ > 1 && cutSide_[i] <= 0) {
      printf("--- %s: Cut %s has function value %.4e with %.3e loss from below (%.1fs)\n",
             __func__, variables_[i].Data(), frac, 1.-npassed/initSignal_, timer->RealTime());
      timer->Reset(); timer->Start();
    }

    //compare function value to previous values, save if best so far
    if(frac > fracBest) {
      fracBest = frac;
      nPassS  = npassed;
      nPassB  = npassedB;
      varBest = i;
      maxOrMin = -1;
      fracStep = (1.-npassed/initSignal_); // total fraction lost after this step
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
    npassed = (cutSide_[i] >= 0) ? GetNumPass(cutSig, sigTreeCut_,sigWeight_) : -1.;

    cutBkg = cutS.Data();
    if(bkgID_ != "") cutBkg+="&&";
    cutBkg += bkgID_.Data();
    npassedB = GetNumPass(cutBkg, bkgTreeCut_,bkgWeight_);
    frac = (cutSide_[i] >= 0) ? EvaluateFunction(npassed, npassedB, function_) : -1.;
    if(verbose_ > 1 && cutSide_[i] >= 0) {
      printf("--- %s: Cut %s has function value %.4e with %.3e loss from above (%.1fs)\n",
             __func__, variables_[i].Data(), frac, 1.-npassed/initSignal_, timer->RealTime());
      timer->Reset();
    }

    //compare function value to previous value, save if best so far
    if(frac > fracBest) {
      fracBest = frac;
      nPassS  = npassed;
      nPassB  = npassedB;
      varBest = i;
      maxOrMin = 1;
      fracStep = 1.-npassed/initSignal_; // total fraction lost after this step
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

  //build the new cut string
  BuildCutString();

  //report the progress if asked to, save the values so far
  if(verbose_ > 0) printf("--- %s: %s = %.4e (signal = %.3e, background = %.3e)\n", __func__, function_.Data(), fracBest, nPassS, nPassB);
  //add efficiencies to the ROC curve
  sigEff_.push_back(nPassS/initSignal_);
  bkgEff_.push_back(1.-nPassB/initBackground_);
  if(verbose_ > 0) printf("--- %s: Signal Efficiency = %.4e (%.3e / %.3e) Background rejection = %.4e (1 - %.3e / %.3e)\n",
                          __func__, nPassS/initSignal_, nPassS, initSignal_, 1.-nPassB/initBackground_, nPassB, initBackground_ );
  cutPath_.push_back(cuts_);
  functionVals_.push_back(fracBest); //record function vs cuts
  EvaluateSpectatorFunctions(nPassS, nPassB); //record spectator function values

  UpdateTrees(); //update the event lists for the trees

  return fracStep;
}

double  CutsetTrainer::GetNumPass(TString selection, TTree* tree, TString func) {
  TString cutstring = selection;
  if(func != "") {
    if(cutstring != "") cutstring = Form("%s*(%s)",func.Data(),selection.Data());
    else cutstring = func;
  }
  TCut cut(cutstring);
  TH1F* hist = new TH1F("hist", "hist", 1, -1e20, 1e20);
  tree->Draw(Form("%s>>hist", variables_[0].Data()), cut, "goff");
  double n = hist->Integral();
  if(n < 0.) {
    printf("!!! %s: Number of events passing the cut are less than 0! n = %.3e, func = %s, cut = %s\n",
           __func__, n, func.Data(), selection.Data());
    n = 1.e-10;
  } else {
    n *= lum_;
  }
  delete hist;
  return n;
}

double CutsetTrainer::EvaluateFunction(double nsignal, double nbackground, TString function) {
  double val = -1.;
  if(function == "s/b") {
    if(nbackground > 0.) val = nsignal/nbackground;
    else return 1.e9*nsignal;
  } else if(function == "significance") {
    if(nsignal > 0. || nbackground > 0.) val = nsignal/sqrt(nsignal + nbackground);
    else val = -1.;
  } else if(function == "limitApprox") {
    if(nbackground > 0.) val = nsignal/sqrt(nbackground)/1.64485;
    else val = 1.e9*nsignal;
  } else if(function == "limit") { //round background down, as n <= nbackground
    if(nbackground > 0.) {
      double p = 0.05; //confidence limit goal
      double tolerance = 0.001; //precision to achieve for goal confidence level
      double scale = 1.; //scale signal until achieve tolerance
      int maxAttempts(20), attempts(0);
      while(abs(val - p) > tolerance && attempts < maxAttempts) { //guess scale factors until close to limit goal
        ++attempts;
        val = ROOT::Math::poisson_cdf((int) (nbackground), nbackground + nsignal*scale); //confidence limit at this value
        if(verbose_ > 5) printf("--- %s: Loop for limit: val = %.3e scale = %.3e\n", __func__, val, scale);
        if(abs(val-p) > tolerance) { //only update if still not succeeding
          scale *= (val/p < 4.) ? (1.-p)/(1.-val) : sqrt(nbackground)*1.64485/nsignal; //if far, start with this scale
          if(attempts == maxAttempts) printf("--- %s: Hit searching attempt limit, %i, with val = %.3e and scale = %.3e\n",
                                            __func__,attempts, val, scale);
        }
      }
      val = 1./scale;
    }
    else val = 1.e9*nsignal;
  } else if(function == "FCDiscovery") { //Feldman-cousins based 5 sigma discovery
    val = significances_.FCDiscovery(nsignal, nbackground);
  } else if(function == "MedianUpperLimit") { //Feldman-cousins based upper limit
    val = nsignal/significances_.FCMedianUpperLimit(nbackground, 0.90); //signal / signal needed for 90% CL
  } else if(function == "MeanUpperLimit") { //Feldman-cousins based upper limit
    val = nsignal/significances_.FCUpperLimit(nbackground, 0.90); //signal / signal needed for 90% CL
  } else if(function == "MeanDiscovery") { //mean significance of 5 sigma
    val = significances_.AverageFiveSigma(nsignal, nbackground);
  } else {
    printf("--- %s: Unknown optimizing function = %s!\n", __func__, function.Data());
    val = -1.;
  }
  return val;
}

void CutsetTrainer::EvaluateSpectatorFunctions(double nsignal, double nbackground) {
  for(unsigned index = 0; index < spectatorFunctions_.size(); ++index) {
    double val = EvaluateFunction(nsignal, nbackground, spectatorFunctions_[index]);
    spectatorFunctionVals_[index].push_back(val);
    if(verbose_ > 0) printf("--- %s: Function %s val = %.3e\n", __func__, spectatorFunctions_[index].Data(), val);
  }
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
  if(verbose_ > 1) printf("--- BuildCutString: Current cut string = %s \n",cuts_.Data());
}

int CutsetTrainer::AddVariable(TString name, int cutMaxOrMin) {

  if(!sigTree_) return 1;
  if(!bkgTree_) return 2;
  if(!sigTree_->GetBranch(name.Data())) return 3;
  if(!bkgTree_->GetBranch(name.Data())) return 4;
  variables_.push_back(name);
  cutSide_.push_back(cutMaxOrMin); //0 = either, > 0 = only max, < 0 means only min
  return 0;
}

int CutsetTrainer::GetMaxMin() {
  Float_t mx[variables_.size()];
  Float_t mn[variables_.size()];
  Float_t vars[variables_.size()];
  TTree* ttemp = sigTreeCut_->CopyTree("");
  if(verbose_ > 2) ttemp->Print();
  ttemp->SetName("ttemp");
  for(int i = 0; i < int(variables_.size()); ++i) {
    mx[i] = -1.e9;
    mn[i] =  1.e9;
    vars[i] = 0.;
    int status = ttemp->SetBranchAddress(variables_[i].Data(), &(vars[i]));
    if(status != 0) {
      if(verbose_ > 0) printf("--- GetMaxMin: Failed to find branch address %s!\n", variables_[i].Data());
      return i+1;
    }
  }
  Long64_t n = ttemp->GetEntriesFast();
  if(verbose_ > 2) printf("--- GetMaxMin: Signal tree has %lld entries\n", n);
  for(Long64_t j = 0; j < n; ++j) {
    ttemp->GetEntry(j);
    for(int i = 0; i < int(variables_.size()); ++i) {
      if     (vars[i] > mx[i]) mx[i] = vars[i];
      else if(vars[i] < mn[i]) mn[i] = vars[i];
      if(verbose_ > 2 && j == 0) printf("--- GetMaxMin: %s entry 0 = %.3e\n", variables_[i].Data(), vars[i]);
    }
  }
  for(int i = 0; i < int(variables_.size()); ++i) {
    if(verbose_ > 2) printf("--- GetMaxMin: %s initial range %.3e - %.3e\n", variables_[i].Data(),
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
  TCanvas* c = new TCanvas("ROC_Canvas","ROC Canvas", 1000, 800);
  roc_->Draw();
  roc_->SetLineColor(kBlue);
  roc_->SetLineWidth(3);
  roc_->SetMarkerStyle(20);
  roc_->SetMarkerSize(0.8);
  roc_->SetTitle("Rectangular Cuts ROC;Signal Efficiency;Background Rejection");
  return c;
}

TCanvas* CutsetTrainer::PlotFunction() {
  int n = sigEff_.size();
  if(n == 0) return NULL;
  Double_t x[n],y[n];
  Double_t min_y(1.e9), max_y(-1.e9);
  for(int i = 0; i < n; ++i) {
    x[i] = sigEff_[i];
    y[i] = functionVals_[i];
    min_y = std::min(y[i], min_y);
    max_y = std::max(y[i], max_y);
  }
  func_ = new TGraph(n,x,y);
  func_->SetName("Function_Graph");
  TCanvas* c = new TCanvas("Function_Canvas","Function Canvas", 1000, 800);
  func_->Draw();
  func_->SetLineColor(kBlue);
  func_->SetLineWidth(3);
  func_->SetMarkerStyle(20);
  func_->SetMarkerSize(0.8);
  func_->SetTitle(Form("Optimizing %s;Signal Efficiency;%s", function_.Data(), function_.Data()));
  TLine* line = new TLine(x[0], max_y, x[n-1], max_y);
  line->SetLineWidth(2);
  line->SetLineColor(kBlue);
  line->Draw("same");

  //Add spectator functions as well
  std::vector<TGraph*> gs;
  if(spectatorFunctions_.size() > 0) {
    const int colors[] = {kRed, kMagenta, kOrange};
    TLegend* leg = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg->AddEntry(func_, function_.Data());
    for(unsigned index = 0; index < spectatorFunctions_.size(); ++index) {
      Double_t xs[n], ys[n];
      Double_t min_ys(1.e9), max_ys(-1.e9);
      for(int i = 0; i < n; ++i) {
        xs[i] = sigEff_[i];
        ys[i] = spectatorFunctionVals_[index][i];
        min_ys = std::min(ys[i], min_ys);
        max_ys = std::max(ys[i], max_ys);
      }
      TGraph* g = new TGraph(n,xs,ys);
      g->SetName(Form("g_%s", spectatorFunctions_[index].Data()));
      g->Draw("PL");
      g->SetLineColor((index >= sizeof(colors)/sizeof(*colors)) ? colors[0] : colors[index]);
      g->SetLineWidth(3);
      g->SetMarkerStyle(21 + index);
      g->SetMarkerSize((index == 0) ? 0.8 : 1.0);
      g->SetLineStyle(kDashed);
      gs.push_back(g);
      leg->AddEntry(g, spectatorFunctions_[index].Data());

      TLine* line_s = new TLine(x[0], max_ys, x[n-1], max_ys);
      line_s->SetLineColor(g->GetLineColor());
      line_s->SetLineWidth(2);
      line_s->SetLineStyle(kDashed);
      line_s->Draw("same");
      max_y = std::max(max_y, max_ys);
      min_y = std::min(min_y, min_ys);
    }
    leg->Draw();
    func_->GetYaxis()->SetRangeUser(0.9*min_y, 1.1*max_y);
  }
  return c;
}

TCanvas* CutsetTrainer::PlotBranchingRatio() {
  int n = sigEff_.size();
  if(n == 0) return NULL;
  Double_t x[n],y[n];
  Double_t min_y(1.e9), max_y(-1.e9), xmin(0.);
  for(int i = 0; i < n; ++i) {
    x[i] = sigEff_[i];
    y[i] = brSig_/functionVals_[i];
    xmin = (y[i] < min_y) ? x[i] : xmin;
    min_y = std::min(y[i], min_y);
    max_y = std::max(y[i], max_y);
  }
  ratio_ = new TGraph(n,x,y);
  ratio_->SetName("BranchingRatio_Graph");
  TCanvas* c = new TCanvas("BranchingRatio_Canvas","Branching Ratio Canvas", 1000, 800);
  ratio_->Draw();
  ratio_->SetLineColor(kBlue);
  ratio_->SetLineWidth(3);
  ratio_->SetMarkerStyle(20);
  ratio_->SetMarkerSize(0.8);
  ratio_->SetTitle(Form("Optimized Branching ratio (%s);Signal Efficiency;Branching ratio (%s)",
                        function_.Data(), function_.Data()));
  TLine* line = new TLine(x[0], min_y, x[n-1], min_y);
  line->SetLineWidth(2);
  line->SetLineColor(kBlue);
  line->Draw("same");
  TLine* line2 = new TLine(xmin, 0.8*min_y, xmin, 1.2*min_y);
  line2->SetLineColor(line->GetLineColor());
  line2->SetLineWidth(2);
  line2->SetLineStyle(kDashed);
  line2->Draw("same");

  //Add spectator functions as well
  std::vector<TGraph*> gs;
  if(spectatorFunctions_.size() > 0) {
    const int colors[] = {kRed, kMagenta, kOrange};
    TLegend* leg = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg->AddEntry(func_, function_.Data());
    double xmin_s(0.);
    for(unsigned index = 0; index < spectatorFunctions_.size(); ++index) {
      Double_t xs[n], ys[n];
      Double_t min_ys(1.e9), max_ys(-1.e9);
      for(int i = 0; i < n; ++i) {
        xs[i] = sigEff_[i];
        ys[i] = brSig_/spectatorFunctionVals_[index][i];
        xmin_s = (ys[i] < min_ys) ? xs[i] : xmin_s;
        min_ys = std::min(ys[i], min_ys);
        max_ys = std::max(ys[i], max_ys);
      }
      TGraph* g = new TGraph(n,xs,ys);
      g->SetName(Form("g_%s", spectatorFunctions_[index].Data()));
      g->Draw("PL");
      g->SetLineColor((index >= sizeof(colors)/sizeof(*colors)) ? colors[0] : colors[index]);
      g->SetLineWidth(3);
      g->SetMarkerStyle(21 + index);
      g->SetMarkerSize((index == 0) ? 0.8 : 1.0);
      g->SetLineStyle(kDashed);
      gs.push_back(g);
      leg->AddEntry(g, spectatorFunctions_[index].Data());

      TLine* line_s = new TLine(x[0], min_ys, x[n-1], min_ys);
      line_s->SetLineColor(g->GetLineColor());
      line_s->SetLineWidth(2);
      line_s->SetLineStyle(kDashed);
      line_s->Draw("same");
      TLine* line_s2 = new TLine(xmin_s, 0.8*min_ys, xmin_s, 1.2*min_ys);
      line_s2->SetLineColor(line_s->GetLineColor());
      line_s2->SetLineWidth(2);
      line_s2->SetLineStyle(kDashed);
      line_s2->Draw("same");

      max_y = std::max(max_y, max_ys);
      min_y = std::min(min_y, min_ys);
    }
    leg->Draw();
    ratio_->GetYaxis()->SetRangeUser(0.7*min_y, 1.1*max_y);
  }
  return c;
}

void CutsetTrainer::WriteResults(TString fname) {
  std::ofstream outfile;
  outfile.open(fname.Data());
  outfile << "#Automatically generated CutsetTrainer result file" << std::endl;
  outfile << "#Path/s:" << function_.Data() << "/d:sigEff/d:bkgEff/d:";
  for(unsigned i = 0; i < spectatorFunctions_.size(); ++i) {
    outfile << spectatorFunctions_[i].Data() << "/d:";
  }
  outfile << std::endl << "#br=" << brSig_ << std::endl
          << "#initSignal=" << initSignal_ << std::endl
          << "#initBackground=" << initBackground_ << std::endl;

  for(unsigned i = 0; i < cutPath_.size(); ++i) {
    outfile << cutPath_[i] << ","
            << functionVals_[i] << ","
            << sigEff_[i] << ","
            << bkgEff_[i];
    for(unsigned j = 0; j < spectatorFunctions_.size(); ++j) {
      outfile << "," << spectatorFunctionVals_[j][i];
    }
    outfile << std::endl;
  }
  outfile.close();
}

bool CutsetTrainer::ReadResults(TString fname) {
  std::ifstream infile;
  infile.open(fname.Data(), std::ifstream::in);
  std::string str;
  int lines = 0;
  //reset the training fields
  cutPath_ = {}; sigEff_ = {}; bkgEff_ ={}; function_ = "";
  functionVals_ = {};
  spectatorFunctions_ = {}; spectatorFunctionVals_ = {};

  //read the input file
  while(std::getline(infile, str)){
    ++lines;
    std::istringstream stream_s(str);

    //read header info
    if(str.find("#Path") != std::string::npos) {
      std::string next_info, sub_info;
      std::getline(stream_s, next_info, ':'); //just Path/s
      std::getline(stream_s, next_info, ':'); //define optimization function
      //get optimizer name
      std::istringstream info_s(next_info);
      std::getline(info_s, sub_info, '/');
      function_ = sub_info.c_str();
      if(verbose_ > 0) printf("--- %s: Optimized function = %s\n", __func__, function_.Data());
      std::getline(stream_s, next_info, ':'); //sigEff
      std::getline(stream_s, next_info, ':'); //bkgEff
      //get spectator functions
      while(std::getline(stream_s, next_info, ':')) {
        if(next_info == "") break;
        std::istringstream info_s2(next_info);
        std::getline(info_s2, sub_info, '/');
        if(verbose_ > 0) printf("--- %s: Spectator function %s identified\n", __func__, sub_info.c_str());
        spectatorFunctions_.push_back(sub_info.c_str());
        spectatorFunctionVals_.push_back({});
      }
      continue;
    }

    //read signal branching ratio used
    if(str.find("#br") != std::string::npos) {
      std::string next_info;
      std::getline(stream_s, next_info, '=');
      std::getline(stream_s, next_info, '=');
      brSig_ = std::stod(next_info);
      continue;
    }

    //read signal norm
    if(str.find("#initSignal") != std::string::npos) {
      std::string next_info;
      std::getline(stream_s, next_info, '=');
      std::getline(stream_s, next_info, '=');
      initSignal_ = std::stod(next_info);
      continue;
    }

    //read background norm
    if(str.find("#initBackground") != std::string::npos) {
      std::string next_info;
      std::getline(stream_s, next_info, '=');
      std::getline(stream_s, next_info, '=');
      initSignal_ = std::stod(next_info);
      continue;
    }

    //ignore comments
    if(str.find("#") != std::string::npos) continue;

    //retrieve the info
    std::string val;
    std::getline(stream_s, val, ',');
    cutPath_.push_back(val.c_str());
    std::getline(stream_s, val, ',');
    functionVals_.push_back(std::stod(val));
    std::getline(stream_s, val, ',');
    sigEff_.push_back(std::stod(val));
    std::getline(stream_s, val, ',');
    bkgEff_.push_back(std::stod(val));
    for(unsigned index = 0; index < spectatorFunctions_.size(); ++index) {
      std::getline(stream_s, val, ',');
      spectatorFunctionVals_[index].push_back(std::stod(val));
    }
  }
  return (lines > 0);
}

void CutsetTrainer::UpdateTrees() {
  TString cutSig = cuts_.Data();
  if(sigID_ != "" && cutSig != "") cutSig += "&&";
  cutSig += sigID_;
  TEventList* list_s = new TEventList("siglist");
  sigTreeCut_->Draw(">>siglist", cutSig.Data(), "goff");
  list_s->SetName("siglist_eval");
  if(list_s->GetN() > 0) sigTreeCut_->SetEventList(list_s);

  TString cutBkg = cuts_.Data();
  if(bkgID_ != "" && cutBkg != "") cutBkg += "&&";
  cutBkg += bkgID_;
  TEventList* list_b = new TEventList("bkglist");
  bkgTreeCut_->Draw(">>bkglist", cutBkg.Data(), "goff");
  list_b->SetName("bkglist_eval");
  if(list_b->GetN() > 0) bkgTreeCut_->SetEventList(list_b);
}

void CutsetTrainer::ReEvaluateFunctions() {
  //clear previous results
  functionVals_ = {};
  for(unsigned index = 0; index < spectatorFunctions_.size(); ++index) {
    if(index >= spectatorFunctionVals_.size()) spectatorFunctionVals_.push_back({});
    else  spectatorFunctionVals_[index] = {};
  }
  //evaluate functions for each signal/background
  for(unsigned index = 0; index < sigEff_.size(); ++index) {
    double nsig = sigEff_[index]*initSignal_;
    double nbkg = (1.-bkgEff_[index])*initBackground_;
    functionVals_.push_back(EvaluateFunction(nsig, nbkg, function_));
    EvaluateSpectatorFunctions(nsig,nbkg);
  }
}
