# CutsetTraining
Training rectangular box cuts on TTrees from ZTauTauHistMaker.

The algorithm is a simple gradient descent, where cuts are selected for small, user defined signal loss steps to optimize a pre-defined function. 

This algorithm only increases cut values, so it will run in `O(s*v*(n_bkg + n_sig*log(n_sig)))`, where s = signal loss steps, v = number of variables, n_bkg = number of events in the background tree, and n_sig = number of events in the signal tree. 

This method naturally ignores variables that aren't useful and has a running time only linear in input variables. Predefined optimization functions are of the form typically of the number of signal over the number of background.

Made for use in ROOT, compiled with ROOT c++ compiler.


An example training script is scripts/train_cms_cutset.C
```
root -l scripts/train_cms_cutset.C
```
