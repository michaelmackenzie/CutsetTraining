

namespace {
  TString folder_ = "../../tmva_training/"; //folder for finding trees
  Double_t signal_fraction_ = 0.5; //fraction of events to train with
  Double_t background_fraction_ = 0.5; //fraction of events to train with
  Int_t verbosity_ = 1;
  TString function_ = "significance"; //function to optimize
  Double_t lum_ = 35.922e3; //luminosity for significance/limits
};


CutsetTrainer* train_cutset(const char* tree_name = "background_ztautau_Z0_mutau_8.tree", double endloss = 0.2, double step = 0.0025,
			    double tolerance = 0.0005, vector<int> signal_ids= {33,34,35, 36, 37, 38, 39}) {

  TString name; //name of output file
  TString x = tree_name;
  TObjArray *tx = x.Tokenize(".");

  printf("Getting Trees\n");

  f = TFile::Open(Form("%s%s",folder_.Data(),tree_name), "READ");

  
  TTree* tree;
  int isData = 0;
  if(x.Contains("mock")) {
    tree = (TTree*) f->Get("background_tree_mock_data");
    isData = 1;
  } else
    tree = (TTree*) f->Get("background_tree");


  if(!tree) {
    printf("Trees not found\n");
    f->ls("");
    return NULL;
  }

  name = "training_";
  
  for(int i = 0; i < tx->GetLast(); ++i) name +=   ((TObjString *)(tx->At(i)))->String();
  // int n = ignore.size();
  // if(n > 0) {
  //   name += "_mva_no_";
  //   if(n > 1)
  //     printf("Ignoring categories: ");
  //   else
  //     printf("Ignoring category: ");
  //   for(int i = 0; i < n; ++i) {
  //     printf("%i%s",ignore[i], (i < n-1) ? ", " : "\n");
  //     name += ignore[i];
  //     if( i < n-1) name += "_";
  //   }
  // }
  TString selection = "";
  if(name.Contains("Z0"))
    selection += "z";
  else if(name.Contains("higgs"))
    selection += "h";
  else {
    printf("Unknown selection! Defaulting to Z0!\n");
    selection += "z";
  }
  
  if(name.Contains("mutau"))
    selection += "mutau";
  else if(name.Contains("etau"))
    selection += "etau";
  else if(name.Contains("emu"))
    selection += "emu";
  else {
    printf("Unknown selection! Default to mutau!\n");
    selection += "mutau";
  }
  selection_ = selection;
  printf("Beginning Training %s with selection %s\n",name.Data(), selection.Data());


  TFile* outf = new TFile("train_tmp.root","RECREATE"); //needed to not only create trees in memory

  Long64_t n = tree->GetEntriesFast();
  outf->cd();
  tree->SetName("Dataset");

  //build signal selection
  TString sCut = "(";
  for(int id : signal_ids) {
    if(sCut.Contains("eventcategory")) sCut += Form(" || eventcategory == %i",id);
    else sCut += Form("eventcategory == %i",id);
  }
  sCut += ")";
  if(sCut == "()") {
    printf("No signals identified!\n");
    return NULL;
  }

  //select tree from data
  printf("*** Identifying signal tree\n");
  TTree* tsignal = tree->CopyTree(sCut.Data(),"");
  tsignal->SetName("signal");

  TString bCut = "!" + sCut; //opposite of signal

  printf("*** Identifying background tree\n");
  TTree* tbackground = tree->CopyTree(bCut.Data(),"");//,n*background_fraction_);
  tbackground->SetName("background");

  CutsetTrainer* trainer = new CutsetTrainer();
  trainer->SetSigTree(tsignal);
  
  trainer->SetBkgTree(tbackground);

  trainer->SetVerbose(verbosity_);

  // trainer->bkgWeight_ = ""; //unit weight
  // trainer->sigWeight_ = "";
  trainer->bkgWeight_ = "fulleventweight"; //unit weight
  trainer->sigWeight_ = "fulleventweight";

  trainer->function_ = function_;
  trainer->lum_ = lum_;
  
  int varStatus = trainer->AddVariable("lepm");
  varStatus = trainer->AddVariable("mtone") || varStatus;
  varStatus = trainer->AddVariable("mttwo") || varStatus;
  varStatus = trainer->AddVariable("leponept") || varStatus;
  varStatus = trainer->AddVariable("leptwopt") || varStatus;
  varStatus = trainer->AddVariable("leppt") || varStatus;
  // varStatus = trainer->AddVariable("njets") || varStatus; //makes very large steps sometimes due to integer values
  varStatus = trainer->AddVariable("lepmestimate") || varStatus;
  varStatus = trainer->AddVariable("onemetdeltaphi") || varStatus;
  varStatus = trainer->AddVariable("leponedeltaphi") || varStatus;
  varStatus = trainer->AddVariable("leptwodeltaphi") || varStatus;
  varStatus = trainer->AddVariable("ht") || varStatus;

  trainer->tolerance_ = tolerance;
  trainer->stepSize_ = step;
  trainer->endLoss_ = endloss;

  printf("signal has %lli entries and background has %lli entries\n",tsignal->GetEntries(),tbackground->GetEntries());
  if(varStatus == 0) {
    int status = trainer->Train();
    printf("trainer status = %i\n", status);
    trainer->BuildCutString();
    printf("Final trainer cuts = %s\n", trainer->cuts_.Data());
    TFile* outroot = new TFile((name+".root").Data(),"RECREATE"); //needed to not only create trees in memory
    TCanvas* croc = trainer->PlotROC();
    croc->Write();
    TCanvas* cfunc = trainer->PlotFunction();
    cfunc->Write();
    outroot->Write();
    ofstream outfile;
    outfile.open ((name+".txt").Data());
    vector<TString> cutPath = trainer->cutPath_;
    for(int i = 0; i < cutPath.size(); ++i) {
      outfile << cutPath[i] << std::endl;
    }
    outfile.close();
  } else 
    printf("Variable selection failed with status = %i\n", varStatus);
  return trainer;
}
