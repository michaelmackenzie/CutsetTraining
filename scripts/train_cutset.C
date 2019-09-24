

namespace {
  const char* fn_par_em = "/mu2e/app/users/mmackenz/mdc2018/MMAnalysis/tmva_training/datasets/e21s721z.tmva_training_trkpatrec.root";
  const char* fn_dar_em = "/mu2e/app/users/mmackenz/mdc2018/MMAnalysis/tmva_training/datasets/e21s721z.tmva_training_calpatrec.root";
  const char* fn_par_ep = "/mu2e/app/users/mmackenz/mdc2018/MMAnalysis/tmva_training/datasets/e41s721z.tmva_training_trkpatrec.root";
  const char* fn_dar_ep = "/mu2e/app/users/mmackenz/mdc2018/MMAnalysis/tmva_training/datasets/e41s721z.tmva_training_calpatrec.root";
};


//-----------------------------------------------------------------------------
// Algorithm    : "trkpatrec" or "calpatrec"
// TrainingMode : "chi2d" or "logfcons"
// BkgWeight    : 0,1,2,3,4  (+100 if use  Z)
//-----------------------------------------------------------------------------
CutsetTrainer* train_cutset(const char* Algorithm = "DAR", const char* Signal = "ep", int bkgWeight = 0, double endloss = 0.2, double step = 0.0025, double tolerance = 0.0005, int oneVar = 0) {

  //  gDirectory->CurrentDirectory()->pwd();
  TString tmvaName, alg, sig;

  alg = Algorithm;
  alg.ToLower();
  sig = Signal;
  sig.ToLower();

  TString outName = "cutset_training_";
  outName += alg;
  outName += "_";
  outName += sig;
  outName += "_BkgWeight";
  outName += bkgWeight;
  TFile* outf = new TFile("train_tmp.root","RECREATE"); //needed to not only create trees in memory

  TString fname;
  if(sig == "em") {
    if      (alg == "par") fname = fn_par_em;
    else if (alg == "dar") fname = fn_dar_em;
  }
  else if (sig == "ep") {
    if      (alg == "par") fname = fn_par_ep;
    else if (alg == "dar") fname = fn_dar_ep;
  } else {
    printf("Combination %s + %s not yet defined, exiting\n", alg.Data(), sig.Data());
    return NULL;
  }

  TFile*  f = TFile::Open(fname.Data(),"READ");

  printf("*** Loading original tree\n");
  f->cd();
  TTree* t = (TTree*) f->Get("tmva_training_tree")->Clone("Dataset");
  Long64_t n = t->GetEntriesFast();
  outf->cd();
  if(!t) {
    printf("Tree not found, exiting\n");
    return NULL;
  }
  t->SetName("Dataset");


  TString sCut = "(p-pmc)<0.25&&(p-pmc)>-0.25";
  if(sig == "ep") sCut += "&&(pmc>87)";
  if(sig == "em") sCut += "&&(pmc>100)&&(tdip>0.57)&&(tdip<1.0)";
  printf("*** Identifying signal tree\n");
  TTree* tsignal = t->CopyTree(sCut.Data(),"",n/12.);//n/12.);
  tsignal->SetName("signal");
  TString bCut = "(p-pmc)>0.7";
  if(sig == "ep") bCut += "&&(pmc>87)";
  if(sig == "em") bCut += "&&(pmc>100)&&(tdip>0.57)&&(tdip<1.0)";
  printf("*** Identifying background tree\n");
  TTree* tbackground = t->CopyTree(bCut.Data(),"",n/2.);
  tbackground->SetName("background");

  CutsetTrainer* trainer = new CutsetTrainer();
   trainer->SetSigTree(tsignal);
   // trainer->SetSigID(sCut);

   trainer->SetBkgTree(tbackground);
   // trainer->SetBkgID(sCut);

  trainer->SetVerbose(1);

  if(bkgWeight == 2) trainer->bkgWeight_ = "max(1.0,exp(2.0*min(p-pmc,3.0)))"; //exp weight
  else               trainer->bkgWeight_ = ""; //unit weight

  int varStatus = trainer->AddVariable("momerr");
  varStatus = trainer->AddVariable("nafract") || varStatus;

  if(oneVar == 0) {
    varStatus = trainer->AddVariable("nactive") || varStatus;
    varStatus = trainer->AddVariable("chi2d") || varStatus;

    varStatus = trainer->AddVariable("t0err") || varStatus;
    varStatus = trainer->AddVariable("d0") || varStatus;
    varStatus = trainer->AddVariable("rmax") || varStatus;
    varStatus = trainer->AddVariable("nda_o_na") || varStatus;
    varStatus = trainer->AddVariable("nza_o_na") || varStatus;
    varStatus = trainer->AddVariable("nma_o_na") || varStatus;
    // varStatus = trainer->AddVariable("tdip") || varStatus;
  }
  trainer->tolerance_ = tolerance;
  trainer->stepSize_ = step;
  trainer->endLoss_ = endloss;
  trainer->sigWeight_="";

  printf("signal has %lli entries and background has %lli entries\n",tsignal->GetEntries(),tbackground->GetEntries());
  if(varStatus == 0) {
    int status = trainer->Train();
    printf("trainer status = %i\n", status);
    trainer->BuildCutString();
    printf("Final trainer cuts = %s\n", trainer->cuts_.Data());
    TFile* outroot = new TFile((outName+".root").Data(),"RECREATE"); //needed to not only create trees in memory
    TCanvas* c = trainer->PlotROC();
    c->Write();
    outroot->Write();
    //    outfile->Close();
    ofstream outfile;
    outfile.open ((outName+".txt").Data());
    vector<TString> cutPath = trainer->cutPath_;
    for(int i = 0; i < cutPath.size(); ++i) {
      outfile << cutPath[i] << std::endl;
    }
    outfile.close();
  } else 
    printf("Variable selection failed with status = %i\n", varStatus);
  return trainer;
}
