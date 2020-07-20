# CutsetTraining
Training rectangular box cuts on TTrees.

The algorithm is a simple gradient descent, where cuts are selected for small, user defined signal loss steps to optimize a pre-defined function. 

This algorithm only increases cut values, so it will run in `O(s*v*(n_bkg + n_sig*log(n_sig)))`, where s = signal loss steps, v = number of variables, n_bkg = number of events in the background tree, and n_sig = number of events in the signal tree. 

This method naturally ignores variables that aren't useful and has a running time only linear in input variables. Predefined optimization functions are of the form typically of the number of signal over the number of background.

Made for use in ROOT, compiled with scons in Mu2e/Offline.

An example training script is scripts/train_cutset.C.
This assumes a TTree from the TRMCMDCAnaModule for various Mu2e tracks.
Edit the file locations in the anonymous namespace of the script for the relevant TTrees.
```
$> root -l 
root> .L scripts/train_cutset.C
root> folder_ = [path/to/tree/dir/]
root> function_ = "[name of optimizing function]"
root> signal_fraction_ = [fraction of events for training]
root> background_fraction_ = [fraction of events for training]
root> lum_ = [luminosity in pb^-1 to assume]
root> CutsetTrainer* cs = train_cutset("[file name]", [fraction of lost signal to end at], \
                                       [signal loss step size], [tolerance on the step size],\
                                       {[list of signal IDs]})
```
