CutsetTrainer* trainer_;
void read_training(TString file = "cutPaths/cutset_training_electron_BkgWeight_0_trainingmode_2_FCUpperLimit_111111.txt") {
  trainer_ = new CutsetTrainer();
  trainer_->SetVerbose(2);
  trainer_->ReadResults(file);
  trainer_->PlotROC();
  trainer_->PlotFunction();
  trainer_->PlotBranchingRatio();
  if(trainer_->initSignal_ < 1e-10){
    trainer_->initSignal_ = 7.411;
    trainer_->initBackground_ = 17.23;
  }
}
