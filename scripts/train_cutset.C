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
  if(test_training_ == 0) {
    cards.push_back(DataCard_t(true ,path_+get_file_name("rpce", 1), "Ana/su2020_RPCAna/Hist"   , "RPC(External)", get_normalization("rpce" , 1,  true), false, true , kGreen-2 ,-10));
    cards.push_back(DataCard_t(false,path_+get_file_name("rpce", 1), "Ana/su2020_RPCAna/Hist"   , "RPC(External)", get_normalization("rpce" , 1,  true), false, true , kGreen-2 ,-9));
    cards.push_back(DataCard_t(true ,path_+get_file_name("rpci", 1), "Ana/su2020_RPCAna/Hist"   , "RPC(Internal)", get_normalization("rpci" , 1,  true), false, true , kGreen-4 ,-10));
    cards.push_back(DataCard_t(false,path_+get_file_name("rpci", 1), "Ana/su2020_RPCAna/Hist"   , "RPC(Internal)", get_normalization("rpci" , 1,  true), false, true , kGreen-4 ,-9));
    cards.push_back(DataCard_t(true ,path_+get_file_name("cosm", 0), "Ana/su2020_CosmicAna/Hist", "Cosmic(lo)"   , get_normalization("cosm" , 1,  true), false, false, kYellow+1, 2));
    cards.push_back(DataCard_t(false,path_+get_file_name("cosm", 0), "Ana/su2020_CosmicAna/Hist", "Cosmic(lo)"   , get_normalization("cosm" , 2,  true), false, false, kYellow+1, 2));
    cards.push_back(DataCard_t(true ,path_+get_file_name("cry3", 0), "Ana/su2020_CosmicAna/Hist", "Cosmic(hi)"   , get_normalization("cry3" , 1,  true), false, false, kOrange+1, 2));
    cards.push_back(DataCard_t(false,path_+get_file_name("cry3", 0), "Ana/su2020_CosmicAna/Hist", "Cosmic(hi)"   , get_normalization("cry3" , 2,  true), false, false, kOrange+1, 2));
    if(doPositron_) {
      cards.push_back(DataCard_t(true ,path_+get_file_name("rmce", 1), "Ana/su2020_RMCAna/Hist"   , "RMC(External)", get_normalization("rmce", 1, false), false, true , kRed+2      ));
      cards.push_back(DataCard_t(false,path_+get_file_name("rmce", 1), "Ana/su2020_RMCAna/Hist"   , "RMC(External)", get_normalization("rmce", 1, false), false, true , kRed+2   , 1));
      cards.push_back(DataCard_t(true ,path_+get_file_name("rmci", 1), "Ana/su2020_RMCAna/Hist"   , "RMC(Internal)", get_normalization("rmci", 1, false), false, true , kRed        ));
      cards.push_back(DataCard_t(false,path_+get_file_name("rmci", 1), "Ana/su2020_RMCAna/Hist"   , "RMC(Internal)", get_normalization("rmci", 1, false), false, true , kRed     , 1));
    }
  }

  if(test_training_ == 2) {
    if(!doPositron_) {
      cards.push_back(DataCard_t(true ,path_+get_file_name("dio_100k" , 1), "Ana/su2020_TrackAna/Hist" , "DIO"     , 3746065./1.e5*get_normalization("dio"  , 1, false), false, true , kViolet-2,0));
    }
    cards.push_back(DataCard_t(true ,path_+get_file_name("cosm_100k", 0), "Ana/su2020_CosmicAna/Hist", "Cosmic(lo)", 429535./1.e5*get_normalization("cosm" , 1,  true), false, false, kYellow+1, 0));
    cards.push_back(DataCard_t(true ,path_+get_file_name("cry3_100k", 0), "Ana/su2020_CosmicAna/Hist", "Cosmic(hi)", 1174968./1.e5*get_normalization("cry3" , 1,  true), false, false, kOrange+1, 0));
  }
  if(doPositron_) {
    cards.push_back(DataCard_t(true ,path_+get_file_name("cpos", 1), "Ana/su2020_ConvAna/Hist" , "#mu^{-}#rightarrow e^{+}", get_normalization("cpos" , 1), true , true , kBlue));
    cards.push_back(DataCard_t(false,path_+get_file_name("cpos", 1), "Ana/su2020_ConvAna/Hist" , "#mu^{-}#rightarrow e^{+}", get_normalization("cpos" , 1), true , true , kBlue, 1));
  } else {
    if(test_training_ < 2) {
      cards.push_back(DataCard_t(true ,path_+get_file_name("cele", 1), "Ana/su2020_ConvAna/Hist" , "#mu^{-}#rightarrow e^{-}", get_normalization("cele" , 1), true , true , kBlue));
      cards.push_back(DataCard_t(false,path_+get_file_name("cele", 2), "Ana/su2020_ConvAna/Hist" , "#mu^{-}#rightarrow e^{-}", get_normalization("cele" , 2), true , true , kBlue));
    } else {
      cards.push_back(DataCard_t(true ,path_+get_file_name("cele_100k", 1), "Ana/su2020_ConvAna/Hist" , "#mu^{-}#rightarrow e^{-}", 378233./1.e5*get_normalization("cele" , 1), true , true , kBlue));
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
  dataplotter_->livetime_[0] = params.one_batch_sec; //dataplotter_->lumi_[0]/pot_1batch*cycle_1batch;
  dataplotter_->livetime_[1] = params.two_batch_sec; //dataplotter_->lumi_[1]/pot_2batch*cycle_2batch;

  return dataplotter_->AddFiles(cards);
}

//add a branch with the normalization included in the weight
int add_normalized_weight(TTree *&tree, float norm) {
  tmp_file_->cd();
  tmp_file_->Add(tree);
  Long64_t nentries = tree->GetEntriesFast();
  Float_t fullweight, eventweight;
  if(!tree->GetBranch("weight")) return 1;
  auto branch = tree->Branch("fullweight", &fullweight, 'F');
  tree->SetBranchAddress("weight", &eventweight);
  for(Long64_t entry = 0; entry < nentries; ++entry) {
    tree->GetEntry(entry,1);
    fullweight = norm*eventweight;
    branch->Fill();
  }
  tree->Write();
  tree->SetBranchStatus("*", 1);
  tmp_file_->Write();
  return 0;
}

//get the signal and background trees, merging them
int initialize_trees(TTree *&signal_tree, TTree *&background_tree) {
  printf("*** Loading trees\n");
  TList* signal_list = new TList;
  TList* background_list = new TList;
  tmp_file_->cd();
  int trk_set = (doPositron_) ? 4024 : 2024;
  TTree* trees[dataplotter_->files_.size()];
  for(unsigned idata = 0; idata < dataplotter_->files_.size(); ++idata) {
    if(do_batch_mode_ > 0 && do_batch_mode_ != 2-dataplotter_->isOneBatch_[idata]) continue;
    if(verbose_ > 0) cout << "*** Loading tree for dataset " << idata << " = " << dataplotter_->labels_[idata].Data()
                          << endl;
    int i_trk_set = trk_set+dataplotter_->setOffsets_[idata];
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
    if(dataplotter_->isBeam_[idata]) norm *= dataplotter_->lumi_[1-dataplotter_->isOneBatch_[idata]];
    else                             norm *= dataplotter_->livetime_[1-dataplotter_->isOneBatch_[idata]];
    if(verbose_ > 0) cout << "*** Adding full event weights to the tree with " << tree_clone->GetEntriesFast()
                          << " using norm = " << norm << endl;
    if(add_normalized_weight(tree_clone, norm)) {
      cout << "!!! Tree for dataset " << idata << " = " << dataplotter_->labels_[idata].Data()
           << " failed to add full event weight!\n";
      return idata+1;
    }
    tmp_file_->Flush();
    tree_clone = (TTree*) tmp_file_->Get(tree_clone->GetName());
    tmp_file_->Add(tree_clone);

    if(dataplotter_->isSignal_[idata]) signal_list    ->Add(tree_clone);
    else                               background_list->Add(tree_clone);
    trees[idata] = tree_clone;
  }
  if(signal_list->GetEntries() > 0) {
    if(verbose_ > 0) cout << "*** Merging signal trees...\n";
    signal_tree = TTree::MergeTrees(signal_list);
    if(!signal_tree) {
      cout << "!!! Signal tree merging failed!\n";
      return 1;
    }
    signal_tree->SetName("signal_tree");
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
    background_tree = TTree::MergeTrees(background_list);
    if(!background_tree) {
      cout << "!!! Background tree merging failed!\n";
      return 1;
    }
    background_tree->SetName("background_tree");
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
                            int bkgWeight = 0, double endloss = 0.2,
                            double step = 0.02, double tolerance = 0.01,
                            TString tag = "111111") {

  test_training_ = 2;
  do_batch_mode_ = 1;
  doPositron_ = positron;

  //initialize the DataPlotter to store files and normalizations
  if(init_dataplotter()) return NULL;

  TString sig = (doPositron_) ? "positron_" : "electron_";
  TString outName = "cutset_training_" + sig + "_BkgWeight_";
  outName += bkgWeight;
  outName += "_" + tag;

  tmp_file_ = new TFile("TMP_train.root","RECREATE"); //needed to not only create trees in memory

  TTree *t_signal, *t_background;
  if(initialize_trees(t_signal, t_background)) return NULL;

  cout << "*** Loaded signal tree with " << t_signal->GetEntriesFast()
       << " entries and background tree with " << t_background->GetEntriesFast()
       << endl;


  CutsetTrainer* trainer = new CutsetTrainer();
  trainer_ = trainer;
  // trainer->function_ = "limitApprox"; //~signal/sqrt(background)
  trainer->function_ = "FCDiscovery";
  trainer->SetSigTree(t_signal);
  // trainer->SetSigID(sCut);

  trainer->SetBkgTree(t_background);
  // trainer->SetBkgID(sCut);

  trainer->SetVerbose(10);

  trainer->bkgWeight_ = "fullweight";
  if(bkgWeight == 2) trainer->bkgWeight_ += "*max(1.0,exp(2.0*min(p-pfront,3.0)))"; //exp weight

  trainer->sigWeight_ = "fullweight";

  //add momentum window of interest
  trainer->sigID_ = (doPositron_) ? "trkqual>-1&&t0>600&&p>89" : "p>102";
  trainer->bkgID_ = (doPositron_) ? "trkqual>-1&&t0>600&&p>89" : "p>102";

  int varStatus = 0;
  varStatus += trainer->AddVariable("p");
  varStatus += trainer->AddVariable("trkqual");
  varStatus += trainer->AddVariable("t0");
  varStatus += trainer->AddVariable("d0");
  varStatus += trainer->AddVariable("tdip");
  varStatus += trainer->AddVariable("pid");

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
  TFile* outroot = new TFile((outName+".root").Data(),"RECREATE"); //needed to not only create trees in memory
  c->Write();
  TCanvas* c2 = trainer->PlotFunction();
  c2->Write(); //also save the function plot
  outroot->Close();
  delete outroot;

  //Write out the selection cuts along with info for that cut
  ofstream outfile;
  outfile.open ((outName+".txt").Data());
  vector<TString> cutPath = trainer->cutPath_;
  for(int i = 0; i < cutPath.size(); ++i) {
    outfile << cutPath[i] << ", "
            << trainer->functionVals_[i] << ", "
            << trainer->sigEff_[i] << ", "
            << trainer->bkgEff_[i]
            << std::endl;
  }
  outfile.close();

  //Remove from disk the temporary file
  gSystem->Exec("rm TMP_train.root;");


  return trainer;
}
