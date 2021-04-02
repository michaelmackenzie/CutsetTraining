#include "su2020/ana/scripts/DataPlotter.C"
#include "su2020/ana/scripts/Parameters.C"
#include "su2020/ana/scripts/init_datasets.C"

TFile* tmp_file_;
CutsetTrainer* trainer_;

DataPlotter* dataplotter_ = 0;
TString path_ = "/mu2e/data/users/mmackenz/su2020/su2020/histograms/"; //path to histogram directory
int test_training_ = 0;

//load data cards into the DataPlotter
int init_dataplotter() {
  std::vector<DataCard_t> cards;
  //constructor:            isOneBatch        fname                       fpath                      label                scale                      isSignal isBeam  color  setOffset
  if(test_training_ >= 0) {
    cards.push_back(DataCard_t(true ,path_+get_file_name("cosm", 0), "Ana/su2020_CosmicAna/Hist", "Cosmic(lo)"   , get_normalization("cosm" , 1,  true), false, false, kYellow+1, 0));
    cards.push_back(DataCard_t(false,path_+get_file_name("cosm", 0), "Ana/su2020_CosmicAna/Hist", "Cosmic(lo)"   , get_normalization("cosm" , 2,  true), false, false, kYellow+1, 1));
    cards.push_back(DataCard_t(true ,path_+get_file_name("cry3", 0), "Ana/su2020_CosmicAna/Hist", "Cosmic(hi)"   , get_normalization("cry3" , 1,  true), false, false, kOrange+1, 0));
    cards.push_back(DataCard_t(false,path_+get_file_name("cry3", 0), "Ana/su2020_CosmicAna/Hist", "Cosmic(hi)"   , get_normalization("cry3" , 2,  true), false, false, kOrange+1, 1));
    cards.push_back(DataCard_t(true ,path_+get_file_name("rpce", 1), "Ana/su2020_RPCAna/Hist"   , "RPC(External)", get_normalization("rpce" , 1, false), false, true , kGreen-2 , 0));
    cards.push_back(DataCard_t(false,path_+get_file_name("rpce", 1), "Ana/su2020_RPCAna/Hist"   , "RPC(External)", get_normalization("rpce" , 2, false), false, true , kGreen-2 , 1));
    cards.push_back(DataCard_t(true ,path_+get_file_name("rpci", 1), "Ana/su2020_RPCAna/Hist"   , "RPC(Internal)", get_normalization("rpci" , 1, false), false, true , kGreen-4 , 0));
    cards.push_back(DataCard_t(false,path_+get_file_name("rpci", 1), "Ana/su2020_RPCAna/Hist"   , "RPC(Internal)", get_normalization("rpci" , 2, false), false, true , kGreen-4 , 1));
    if(doPositron_) {
      cards.push_back(DataCard_t(true ,path_+get_file_name("rmce", 1), "Ana/su2020_RMCAna/Hist" , "RMC(External)", get_normalization("rmce" , 1, false), false, true , kRed+2   , 0));
      cards.push_back(DataCard_t(false,path_+get_file_name("rmce", 1), "Ana/su2020_RMCAna/Hist" , "RMC(External)", get_normalization("rmce" , 1, false), false, true , kRed+2   , 1));
      cards.push_back(DataCard_t(true ,path_+get_file_name("rmci", 1), "Ana/su2020_RMCAna/Hist" , "RMC(Internal)", get_normalization("rmci" , 1, false), false, true , kRed     , 0));
      cards.push_back(DataCard_t(false,path_+get_file_name("rmci", 1), "Ana/su2020_RMCAna/Hist" , "RMC(Internal)", get_normalization("rmci" , 1, false), false, true , kRed     , 1));
    }
  }

  if(test_training_ == 2) {
    if(!doPositron_) {
      cards.push_back(DataCard_t(true ,path_+get_file_name("dio_100k" , 1), "Ana/su2020_TrackAna/Hist" , "DIO"     , 3746065./1.e5*get_normalization("dio"  , 1, false), false, true , kViolet-2,0));
      cards.push_back(DataCard_t(false,path_+get_file_name("dio_100k" , 1), "Ana/su2020_TrackAna/Hist" , "DIO"     , 3746065./1.e5*get_normalization("dio"  , 2, false), false, true , kViolet-2,1));
    }
    // cards.push_back(DataCard_t(true ,path_+get_file_name("cosm_100k", 0), "Ana/su2020_CosmicAna/Hist", "Cosmic(lo)", 429535./1.e5*get_normalization("cosm" , 1,  true), false, false, kYellow+1, 0));
    // cards.push_back(DataCard_t(true ,path_+get_file_name("cry3_100k", 0), "Ana/su2020_CosmicAna/Hist", "Cosmic(hi)", 1174968./1.e5*get_normalization("cry3" , 1,  true), false, false, kOrange+1, 0));
  } else if(test_training_ < 2) {
    if(!doPositron_) {
      cards.push_back(DataCard_t(true ,path_+get_file_name("dio" , 1), "Ana/su2020_TrackAna/Hist" , "DIO"     , get_normalization("dio"  , 1, false), false, true , kViolet-2,0));
      cards.push_back(DataCard_t(false,path_+get_file_name("dio" , 1), "Ana/su2020_TrackAna/Hist" , "DIO"     , get_normalization("dio"  , 2, false), false, true , kViolet-2,1));
    }
  }

  if(doPositron_) {
    cards.push_back(DataCard_t(true ,path_+get_file_name("cpos", 1), "Ana/su2020_ConvAna/Hist" , "#mu^{-}#rightarrow e^{+}", get_normalization("cpos" , 1), true , true , kBlue, 0));
    cards.push_back(DataCard_t(false,path_+get_file_name("cpos", 1), "Ana/su2020_ConvAna/Hist" , "#mu^{-}#rightarrow e^{+}", get_normalization("cpos" , 1), true , true , kBlue, 1));
  } else {
    if(test_training_ < 2) {
      cards.push_back(DataCard_t(true ,path_+get_file_name("cele", 1), "Ana/su2020_ConvAna/Hist" , "#mu^{-}#rightarrow e^{-}", get_normalization("cele" , 1), true , true , kBlue, 0));
      cards.push_back(DataCard_t(false,path_+get_file_name("cele", 2), "Ana/su2020_ConvAna/Hist" , "#mu^{-}#rightarrow e^{-}", get_normalization("cele" , 2), true , true , kBlue, 1));
    } else {
      cards.push_back(DataCard_t(true ,path_+get_file_name("cele_100k", 1), "Ana/su2020_ConvAna/Hist" , "#mu^{-}#rightarrow e^{-}", 378233./1.e5*get_normalization("cele" , 1), true , true , kBlue, 0));
      cards.push_back(DataCard_t(false,path_+get_file_name("cele_100k", 1), "Ana/su2020_ConvAna/Hist" , "#mu^{-}#rightarrow e^{-}", 378233./1.e5*get_normalization("cele" , 2), true , true , kBlue, 1));
    }
  }
  if(dataplotter_) delete dataplotter_;
  dataplotter_ = new DataPlotter();
  Parameters params(doPositron_); //Running parameters
  dataplotter_->lumi_[0]     = params.one_batch_pot*(do_batch_mode_ != 2); //one batch
  dataplotter_->lumi_[1]     = params.two_batch_pot*(do_batch_mode_ != 1); //two batch
  // const double cycle_1batch  = 1.33; //seconds
  // const double pot_1batch    = 4.e12; //pot/cycle
  // const double cycle_2batch  = 1.4; //seconds
  // const double pot_2batch    = 8.e12; //pot/cycle
  dataplotter_->livetime_[0] = params.one_batch_sec*(do_batch_mode_ != 2); //dataplotter_->lumi_[0]/pot_1batch*cycle_1batch;
  dataplotter_->livetime_[1] = params.two_batch_sec*(do_batch_mode_ != 1); //dataplotter_->lumi_[1]/pot_2batch*cycle_2batch;

  return dataplotter_->AddFiles(cards);
}

//merge a list of trees
TTree* merge_trees(TList* list) {
  TTree* tree = 0;
  int merge_method = (list->GetEntries() <= 1) ? 1 : 0;
  if(merge_method == 0) {
    tree = TTree::MergeTrees(list);
    tmp_file_->Flush();
    if(tree) tree->SetName("tree_merge");
    else return NULL;
  } else if(merge_method == 1) {
    if(list->GetEntries() == 0) return NULL;
    tree = (TTree*) list->At(0);
    tree = tree->CloneTree();
    tree->SetName("tree_merge");
    for(int itree = 1; itree < list->GetEntries(); ++itree) {
      tree->CopyEntries((TTree*) list->At(itree));
    }
  }
  return tree;
}

//add a branch with the normalization included in the weight
int add_normalized_weight(TTree *&tree, float norm, int idata) {
  tmp_file_->cd();
  tmp_file_->Add(tree);
  auto cachesize = 50000000U; //100MB
  tree->SetCacheSize(cachesize);
  tree->AddBranchToCache("*",true);
  tree->LoadBaskets();
  Long64_t nentries = tree->GetEntriesFast();
  Float_t fullweight, eventweight;
  if(!tree->GetBranch("weight")) return 1;
  auto branch = tree->Branch("fullweight", &fullweight, 'F');
  tree->SetBranchAddress("weight", &eventweight);
  auto branch_cat = tree->Branch("category", &idata, 'I');
  for(Long64_t entry = 0; entry < nentries; ++entry) {
    tree->GetEntry(entry,1);
    fullweight = norm*eventweight;
    branch->Fill();
    branch_cat->Fill(); //to be able to track down the originating file
  }
  tree->SetBranchStatus("*", 1);
  tree->Write();
  // tmp_file_->Write();
  return 0;
}

//get the signal and background trees, merging them
int initialize_trees(TTree *&signal_tree, TTree *&background_tree, int trk_set) {
  printf("*** Loading trees\n");
  TList* signal_list = new TList;
  TList* background_list = new TList;
  tmp_file_->cd();
  vector<unsigned> tree_indices;
  for(unsigned idata = 0; idata < dataplotter_->files_.size(); ++idata) {
    if(do_batch_mode_ > 0 && do_batch_mode_ != 2-dataplotter_->isOneBatch_[idata]) continue;
    if(verbose_ > 0) cout << "*** Loading tree for dataset " << idata << " = " << dataplotter_->labels_[idata].Data()
                          << endl;
    int i_trk_set = trk_set + dataplotter_->setOffsets_[idata];
    TTree* tree = (TTree*) dataplotter_->files_[idata]->Get(Form("trk_%i/tracktree_%i", i_trk_set, i_trk_set));
    if(!tree) {
      cout << "!!! Tree for dataset " << idata << " = " << dataplotter_->labels_[idata].Data()
           << " is not found!\n";
      return idata+1;
    }
    tmp_file_->cd();
    TTree* tree_clone = tree->CloneTree(); tree_clone->SetName(Form("tree_%i", idata));
    tmp_file_->Add(tree_clone);

    float norm = dataplotter_->scales_[idata];
    if(dataplotter_->isBeam_[idata]) norm *= dataplotter_->lumi_    [1-dataplotter_->isOneBatch_[idata]];
    else                             norm *= dataplotter_->livetime_[1-dataplotter_->isOneBatch_[idata]];
    if(verbose_ > 0) cout << "*** Adding full event weights to the tree with " << tree_clone->GetEntriesFast()
                          << " using norm = " << norm << endl;
    if(add_normalized_weight(tree_clone, norm, idata)) {
      cout << "!!! Tree for dataset " << idata << " = " << dataplotter_->labels_[idata].Data()
           << " failed to add full event weight!\n";
      return idata+1;
    }
    tree_indices.push_back(idata);
  }

  tmp_file_->Flush();

  //make the list of trees
  for(unsigned idata : tree_indices) {
    TTree* tree_clone = (TTree*) tmp_file_->Get(Form("tree_%i", idata));
    tmp_file_->Add(tree_clone);
    if(dataplotter_->isSignal_[idata]) signal_list    ->Add(tree_clone);
    else                               background_list->Add(tree_clone);
  }

  //merge the list of trees
  if(signal_list->GetEntries() > 0) {
    if(verbose_ > 0) cout << "*** Merging signal trees...\n";
    signal_tree = merge_trees(signal_list);
    if(!signal_tree) {
      cout << "!!! Signal tree merging failed!\n";
      return 1;
    }
    if(verbose_ > 0) cout << "*** Merged signal trees!\n";
    signal_tree->SetName("signal_tree");
    tmp_file_->Add(signal_tree);
    signal_tree->Write();
    if(!signal_tree->GetBranch("fullweight")) {
      cout << "!!! Full event weight branches not properly defined!\n";
      signal_tree->Print();
      return 3;
    } else {
      for(auto obj : *signal_list) delete obj;
    }
  } else return 1;

  if(background_list->GetEntries() > 0) {
    if(verbose_ > 0) cout << "*** Merging background trees...\n";
    background_tree = merge_trees(background_list);
    if(!background_tree) {
      cout << "!!! Background tree merging failed!\n";
      return 1;
    }
    if(verbose_ > 0) cout << "*** Merged background trees!\n";
    background_tree->SetName("background_tree");
    tmp_file_->Add(background_tree);
    background_tree->Write();
    if(!background_tree->GetBranch("fullweight")) {
      cout << "!!! Full event weight branches not properly defined!\n";
      background_tree->Print();
      return 4;
    } else {
      for(auto obj : *background_list) delete obj;
    }
  } else return 2;
  return 0;
}

CutsetTrainer* train_cutset(bool positron = false,
                            int bkgWeight = 0, double endloss = 0.8,
                            double step = 0.01, double tolerance = 0.001,
                            TString tag = "111111",
                            TString function = "MeanUpperLimit"/*"MeanDiscovery" "MeanUpperLimit" "MedianUpperLimit" "FCDiscovery"*/) {

  test_training_ = 2;
  do_batch_mode_ = 0;
  doPositron_ = positron;

  //initialize the DataPlotter to store files and normalizations
  if(init_dataplotter()) return NULL;

  TString sig = (doPositron_) ? "positron" : "electron";
  TString outName = "cutset_training_" + sig + "_BkgWeight_";
  outName += bkgWeight;
  outName += "_trainingmode_";
  outName += test_training_;
  outName += "_" + function;
  if(do_batch_mode_ > 0) {
    outName += "_batchmode_";
    outName += do_batch_mode_;
  }
  if(tag != "") {
    outName += "_" + tag;
  }

  tmp_file_ = new TFile("TMP_train.root","RECREATE", "", 101); //needed to not only create trees in memory

  TTree *t_signal, *t_background;
  int trk_set = (doPositron_) ? 4024 : 2024;
  if(initialize_trees(t_signal, t_background, trk_set)) return NULL;

  cout << "*** Loaded signal tree with " << t_signal->GetEntriesFast()
       << " entries and background tree with " << t_background->GetEntriesFast()
       << endl;


  CutsetTrainer* trainer = new CutsetTrainer();
  trainer_ = trainer;
  // trainer->function_ = "limitApprox"; //~signal/sqrt(background)
  trainer->function_ = function; //training function
  std::vector<TString> spectators;
  if(function == "FCDiscovery")        spectators = {"FCUpperLimit", "MeanDiscovery"                 };
  else if(function == "FCUpperLimit")  spectators = {"FCDiscovery" , "MeanDiscovery"                 };
  else if(function == "MeanDiscovery") spectators = {"FCDiscovery" , "FCUpperLimit"                  };
  else                                 spectators = {"FCDiscovery" , "FCUpperLimit" , "MeanDiscovery"};
  trainer->SetSpectatorFunctions(spectators);

  trainer->brSig_ = br_conv_; //for plotting purposes

  trainer->SetSigTree(t_signal);
  // trainer->SetSigID(sCut);

  trainer->SetBkgTree(t_background);
  // trainer->SetBkgID(sCut);

  trainer->SetVerbose(1);

  trainer->bkgWeight_ = "fullweight";
  if(bkgWeight == 2) trainer->bkgWeight_ += "*max(1.0,exp(2.0*min(p-pfront,3.0)))"; //exp weight

  trainer->sigWeight_ = "fullweight";

  //add momentum window of interest
  trainer->sigID_ = (doPositron_) ? "p>89" : "trkqual>0.02&&p>102.5&&weight<1.e6";
  trainer->bkgID_ = (doPositron_) ? "p>89" : "trkqual>0.02&&p>102.5&&weight<1.e6";

  int varStatus = 0;
  varStatus += trainer->AddVariable("p");
  varStatus += trainer->AddVariable("trkqual", -1); //only cut from below
  varStatus += trainer->AddVariable("t0");
  varStatus += trainer->AddVariable("d0");
  varStatus += trainer->AddVariable("tdip");
  varStatus += trainer->AddVariable("pid", -1); //only cut from below

  // varStatus += trainer->AddVariable("momerr");
  // varStatus += trainer->AddVariable("nafract");
  // varStatus += trainer->AddVariable("nactive");
  // varStatus += trainer->AddVariable("chi2d");
  // varStatus += trainer->AddVariable("t0err");
  // varStatus += trainer->AddVariable("rmax");
  // varStatus += trainer->AddVariable("nda_o_na");
  // varStatus += trainer->AddVariable("nza_o_na");
  // varStatus += trainer->AddVariable("nma_o_na");

  trainer->tolerance_ = tolerance;
  trainer->stepSize_ = step;
  trainer->endLoss_ = endloss;

  if(varStatus) {
    cout << "!!! Trainer variable definition returned a total status of " << varStatus
         << " --> exiting!\n";
    return trainer;
  }

  cout << "*** Trainer variables defined, beginning training!\n";
  int status = trainer->Train();

  //Print final cuts
  printf("trainer status = %i\n", status);
  trainer->BuildCutString();
  printf("Final trainer cuts = %s\n", trainer->cuts_.Data());

  //Write out the ROC plot
  TCanvas* c = trainer->PlotROC();
  TFile* outroot = new TFile((outName+".root").Data(),"RECREATE");
  c->Write(); c->Print(Form("figures/%s_ROC.png", outName.Data()));
  TCanvas* c2 = trainer->PlotFunction();
  //also save the function plot
  c2->Write(); c2->Print(Form("figures/%s_Func.png", outName.Data()));
  TCanvas* c3 = trainer->PlotBranchingRatio();
  //also save the branching ratio plot
  c3->Write(); c3->Print(Form("figures/%s_Br.png", outName.Data()));
  outroot->Close();
  delete outroot;

  //Write out the selection cuts along with info for that cut
  //FIXME: Add output/input of cutpath info to the trainer
  trainer->WriteResults(outName+".txt");

  //Remove from disk the temporary file
  gSystem->Exec("rm TMP_train.root;");


  return trainer;
}
