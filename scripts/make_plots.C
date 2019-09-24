
namespace {
  const char* fn_dar_roc_em_0_cut = "/mu2e/app/users/mmackenz/mdc2018/MMAnalysis/cutset_training/scripts/rootFiles/cutset_training_dar_em_BkgWeight0.root";
  const char* fn_dar_roc_em_2_cut = "/mu2e/app/users/mmackenz/mdc2018/MMAnalysis/cutset_training/scripts/rootFiles/cutset_training_dar_em_BkgWeight2.root";
  const char* fn_dar_roc_ep_0_cut = "/mu2e/app/users/mmackenz/mdc2018/MMAnalysis/cutset_training/scripts/rootFiles/cutset_training_dar_ep_BkgWeight0.root";
  const char* fn_dar_roc_ep_2_cut = "/mu2e/app/users/mmackenz/mdc2018/MMAnalysis/cutset_training/scripts/rootFiles/cutset_training_dar_ep_BkgWeight2.root";

  const char* fn_dar_roc_em_0_mva = "/mu2e/app/users/mmackenz/mdc2018/MMAnalysis/tmva_training/dar_e-_chi2d_BkgWeight0_tmva.root";
  const char* fn_dar_roc_em_2_mva = "/mu2e/app/users/mmackenz/mdc2018/MMAnalysis/tmva_training/dar_e-_chi2d_BkgWeight2_tmva.root";
  const char* fn_dar_roc_ep_0_mva = "/mu2e/app/users/mmackenz/mdc2018/MMAnalysis/tmva_training/dar_e+_chi2d_BkgWeight0_tmva.root";
  const char* fn_dar_roc_ep_2_mva = "/mu2e/app/users/mmackenz/mdc2018/MMAnalysis/tmva_training/dar_e+_chi2d_BkgWeight2_tmva.root";

  const char* fn_par_em = "/mu2e/app/users/mmackenz/mdc2018/MMAnalysis/tmva_training/datasets/e21s721z.tmva_training_trkpatrec.root";
  const char* fn_dar_em = "/mu2e/app/users/mmackenz/mdc2018/MMAnalysis/tmva_training/datasets/e21s721z.tmva_training_calpatrec.root";
  const char* fn_par_ep = "/mu2e/app/users/mmackenz/mdc2018/MMAnalysis/tmva_training/datasets/e41s721z.tmva_training_trkpatrec.root";
  const char* fn_dar_ep = "/mu2e/app/users/mmackenz/mdc2018/MMAnalysis/tmva_training/datasets/e41s721z.tmva_training_calpatrec.root";

  const char* fn_stn_em = "/mu2e/app/users/mmackenz/mdc2018/MMAnalysis/tmva_training/datasets/e21s721z.tcp_ana.hist";
  const char* fn_stn_ep = "/mu2e/app/users/mmackenz/mdc2018/MMAnalysis/tmva_training/datasets/e41s721z.tcp_ana.hist";
};

//Double-sided Crystal Ball function for use with TF1
Double_t crystal(Double_t *x, Double_t *par) {

  Double_t norm   = par[0];
  Double_t mean   = par[1];
  Double_t sigma  = par[2];
  Double_t alpha1 = par[3];
  Double_t n1     = par[4];
  Double_t alpha2 = par[5];
  Double_t n2     = par[6];

  Double_t res = 0.;
  //require non zero sigma and tails on opposite sides of the gaussian, and decaying tails
  //redefined alphas to both be positive numbers, number of sigmas from mean
  if(sigma >= 0. && (alpha1 > 0.) && (alpha2 > 0.) && n1 > 0. && n2 > 0.) {
    Double_t arg = (x[0]-mean)/sigma;

    if(arg < -1.*alpha1) { //more than alpha1 sigmas below mean
      res = pow(n1/(alpha1),n1)*exp(-alpha1*alpha1/2.)*pow((n1/alpha1-alpha1-arg),-n1);
    } else if(arg > alpha2) { //more than alpha2 sigmas above mean
      res = pow(n2/alpha2,n2)*exp(-alpha2*alpha2/2.)*pow((n2/alpha2-alpha2+arg),-n2);
    } else { //else, gaussian
      res = exp(-arg*arg/2.);
    }
  } else if(x[0]==0){
    printf("Failed test, vals: %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n",
	   x[0], norm, mean, sigma, alpha1, n1, alpha2, n2);
  }
  res *= norm;

  return res;
}


TCanvas* plot_dpf(const char* Signal = "ep", const char* selection = "") {
  TString sig = Signal;
  sig.ToLower();
  TFile *fp, *fd;
  if(sig == "ep") {
    fp = TFile::Open(fn_par_ep,"READ");
    fd = TFile::Open(fn_dar_ep,"READ");
  } else if(sig == "em") {
    fp = TFile::Open(fn_par_em,"READ");
    fd = TFile::Open(fn_dar_em,"READ");
  } else
    return NULL;

  if(!fd || !fp) {
    printf("Files not found, exiting\n");
    return NULL;
  }

  TCanvas* cDpF = new TCanvas("cDpF","cDpF");
  TH1F *hPAR = new TH1F("hPAR", "PAR", 100, -5., 5.);
  TH1F *hDAR = new TH1F("hDAR", "DAR", 100, -5., 5.);

  TTree* tp = (TTree*) fp->Get("tmva_training_tree");
  tp->SetName("PAR");
  printf("Getting PAR tree\n");
  tp->Draw("(p-pmc)>>hPAR", selection);
  printf("Drawing PAR tree\n");
  hPAR = (TH1F*) gDirectory->Get("hPAR");

  TTree* td = (TTree*) fd->Get("tmva_training_tree");
  td->SetName("DAR");

  if(!td || !tp) {
    printf("Trees not found, exiting\n");
    return NULL;
  }

  td->Draw("(p-pmc)>>hDAR", selection);
  
  hDAR = (TH1F*) gDirectory->Get("hDAR");

  gStyle->SetOptStat(1111);

  hDAR->SetLineWidth(2);
  hDAR->SetLineColor(kBlue);
  hDAR->Draw();

  hPAR->SetLineWidth(2);
  hPAR->SetLineColor(kRed);
  hPAR->Draw("sames");
  TPaveStats *st = (TPaveStats*) hPAR->FindObject("stats");
  if(st) {
    st->SetY1NDC(0.5); //new y start position
    st->SetY2NDC(0.7); //new y end position
  }
  //  hDAR->SetAxisRange(-5.,5.);
  hDAR->SetTitle("Track Momentum - MC Momentum");
  hDAR->SetXTitle("Reconstructed P - MC P (MeV)");
  hDAR->SetYTitle(Form("Entries / %.1f MeV",hDAR->GetBinWidth(1)));

  cDpF->SetLogy();
  cDpF->SetGridx();
  cDpF->SetGridy();

  TLegend* legend = new TLegend();
  legend->AddEntry(hDAR,"DAR");
  legend->AddEntry(hPAR,"PAR");
  legend->Draw();
  return cDpF;
}

TCanvas* plot_roc(const char* Signal = "ep", const char* Resolver = "dar", int bkg_weight = 0) {
  TString sig = Signal;
  TString res = Resolver;
  sig.ToLower();
  res.ToLower();
  TFile *fMVA, *fCut;
  if(res != "dar") {
    printf("par not implemented yet\n");
    return NULL;
  }
  if(sig == "em") {
    if(bkg_weight == 0) {
      fMVA = (TFile*) TFile::Open(fn_dar_roc_em_0_mva,"READ")->Get("tmva_dar_e-_chi2d_BkgWeight0_tmva");
      fCut = TFile::Open(fn_dar_roc_em_0_cut,"READ");
    }
    if(bkg_weight == 2) {
      fMVA = (TFile*) TFile::Open(fn_dar_roc_em_2_mva,"READ")->Get("tmva_dar_e-_chi2d_BkgWeight2_tmva");
      fCut = TFile::Open(fn_dar_roc_em_2_cut,"READ");
    }
  } else if(sig == "ep") {
    if(bkg_weight == 0) {
      fMVA = (TFile*) TFile::Open(fn_dar_roc_ep_0_mva,"READ")->Get("tmva_dar_e+_chi2d_BkgWeight0_tmva");
      fCut = TFile::Open(fn_dar_roc_ep_0_cut,"READ");
    }
    if(bkg_weight == 2) {
      fMVA = (TFile*) TFile::Open(fn_dar_roc_ep_2_mva,"READ")->Get("tmva_dar_e+_chi2d_BkgWeight2_tmva");
      fCut = TFile::Open(fn_dar_roc_ep_2_cut,"READ");
    }
  } else
    return NULL;

  if(!fMVA || !fCut) {
    printf("Files not found, exiting\n");
    return NULL;
  }

  TCanvas* c = (TCanvas*) fCut->Get("ROC Canvas");
  c->SetName("ROC_Canvas");
  c->SetWindowSize(1000,600);
  TGraph* g = (TGraph*) c->GetPrimitive("Graph");
  c->Draw();
  TFile* f = (TFile*) ((TFile*) fMVA->Get("Method_MLP"))->Get("MLP");
  TH1D* hMVA = (TH1D*) f->Get("MVA_MLP_rejBvsS");
  hMVA->Draw("same C");
  hMVA->SetLineColor(kRed);
  hMVA->SetLineWidth(2);
  gStyle->SetOptStat(0);
  g->SetTitle("ROC;Signal Efficiency;Background Rejection");
  c->SetGridx();
  c->SetGridy();

  TLegend* legend = new TLegend();
  legend->AddEntry(hMVA,"MLP");
  legend->AddEntry(g,"Cut");
  legend->Draw();
  c->Print(Form("%s_%s_BkgWeight_%i_ROC.png",res.Data(),sig.Data(),bkg_weight));
  return c;
}

TCanvas* plot_var(const char* var = "(p-pmc)", const char* Signal = "ep", const char* units = "MeV", 
		  int logPlot = 1, int normalize = 1,
		  const char* resolver = "dar", const char* selection = "") {
  TString sig = Signal;
  TString res = resolver;
  sig.ToLower();
  res.ToLower();
  TFile *f;
  if(sig == "ep") {
    if(res == "dar")
      f = TFile::Open(fn_par_ep,"READ");
    else if(res == "par")
      f = TFile::Open(fn_dar_ep,"READ");
  } else if(sig == "em") {
    if(res == "dar")
      f = TFile::Open(fn_par_em,"READ");
    else if(res == "par")
      f = TFile::Open(fn_dar_em,"READ");
  } else
    return NULL;

  if(!f) {
    printf("File not found, exiting\n");
    return NULL;
  }

  TCanvas* c = new TCanvas("c","c");
  TH1F *hBackground;// = new TH1F("hBackground", "Background", 100, -5., 5.);
  TH1F *hSignal;// = new TH1F("hSignal", "Signal", 100, -5., 5.);
  TString bkgCut = "(p-pmc)>0.7";
  TString sel = selection;
  if(TString(var) == "chi2d") {
    hSignal = new TH1F("hSignal","hSignal",100,0.,10.);
    hBackground = new TH1F("hBackground","hBackground",100,0.,10.);
  }
  if(TString(var) == "(p-pmc)") {
    hSignal = new TH1F("hSignal","hSignal",100,-5.,5.);
    hBackground = new TH1F("hBackground","hBackground",100,-5.,5.);
  }

  if(sel != "") 
    bkgCut += "&&";
  bkgCut += sel;
  TString sigCut = "(p-pmc)>-0.25&&(p-pmc)<0.25";
  if(sel != "") 
    sigCut += "&&";
  sigCut += sel;
  TString varDraw = var;
  
  TTree* t = (TTree*) f->Get("tmva_training_tree");
  printf("Drawing Background tree\n");
  t->Draw((varDraw+">>hBackground").Data(), bkgCut.Data());
  hBackground = (TH1F*) gDirectory->Get("hBackground");
  t->Draw((varDraw+">>hSignal").Data(), sigCut.Data());
  printf("Drawing Signal tree\n");
  hSignal = (TH1F*) gDirectory->Get("hSignal");

  if(!hSignal || !hBackground) {
    printf("Histograms not found, exiting\n");
    return NULL;
  }

  gStyle->SetOptStat(0);

  if(normalize > 0) {
    hSignal->Scale(1./hSignal->Integral());
    hBackground->Scale(1./hBackground->Integral());
  }
  hSignal->SetLineWidth(2);
  hSignal->SetLineColor(kBlue);
  hSignal->SetFillColor(kAzure-4);
  hSignal->Draw("hist");

  hBackground->SetLineWidth(2);
  hBackground->SetLineColor(kRed);
  hBackground->SetFillStyle(3002);
  hBackground->SetFillColor(kRed);
  hBackground->Draw("hist sames");

  TString varTitle = var;
  hSignal->SetTitle((varTitle+" Distribution").Data());
  if(TString(units) != "") {
    varTitle += " (";
    varTitle += units;
    varTitle += ")";
  }
  hSignal->SetXTitle(varTitle.Data());
  hSignal->SetYTitle(Form("Entries / %.1e %s",hSignal->GetBinWidth(1),units));
  if(hSignal->GetMaximum() < hBackground->GetMaximum())
    hSignal->SetMaximum(hBackground->GetMaximum());
  if(logPlot)  c->SetLogy();
  c->SetGridx();
  c->SetGridy();
  // c->SetFillColor(kGray);
  // c->GetFrame()->SetFillColor(kWhite);
  TLegend* legend = new TLegend();
  legend->AddEntry(hSignal,"Signal");
  legend->AddEntry(hBackground,"Background");
  legend->Draw();
  c->Print(Form("%s_%s_plot.png",var,Signal));
  return c;
}

/*
  dataset      | MLP Cut Val | Efficiency | Cut Number | Efficiency
  dar_e+_bkg_0 |  0.25       | 90.1 %     |     6      |  
  dar_e-_bkg_0 |  0.35       | 90.1 %     |            |
  dar_e+_bkg_2 |  0.70       | 89.9 %     |            |
  dar_e-_bkg_2 |  0.87       | 89.7 %     |            |


 */
TCanvas* compare_mva(const char* Signal = "ep",
		     int logPlot = 1, int fitCrystal = 0,
		     const char* resolver = "dar",
		     const char* cuts = "", int bkg_weight = 0) {

  TString sig = Signal;
  TString res = resolver;
  TString cut = cuts;
  sig.ToLower();
  res.ToLower();
  TFile *f_mva, *f_cut, *f_stn, *f_tot;
  if (bkg_weight != 0 && bkg_weight != 2) {
    printf("Background weight unknown\n");
    return NULL;
  }
  if(sig == "ep") {
    f_stn = TFile::Open(fn_stn_ep,"READ");
    if(res == "dar")
      f_cut = TFile::Open(fn_dar_ep,"READ");
    else
      f_cut = TFile::Open(fn_par_ep,"READ");
  } else if(sig == "em") {
    f_stn = TFile::Open(fn_stn_em,"READ");
    if(res == "dar")
      f_cut = TFile::Open(fn_dar_em,"READ");
    else
      f_cut = TFile::Open(fn_par_em,"READ");
  }
  else {
    printf("Signal not defined\n");
    return NULL;
  }
  int ihist = 130;
  if(res == "dar") ihist = 235;
  if(sig == "ep")  ihist += 2;
  if(bkg_weight == 2) ihist += 1; 
  
  f_mva = (TFile*) ((TFile*) ((TFile*) ((TFile*) f_stn->Get("Ana"))->Get("TrackComp"))->Get("Hist"))->Get(Form("trk_%i",ihist));
  int itot = (res == "dar") ? 207 : 107; //ntrack == 1 to match ttree writing
  f_tot = (TFile*) ((TFile*) ((TFile*) ((TFile*) f_stn->Get("Ana"))->Get("TrackComp"))->Get("Hist"))->Get(Form("trk_%i",itot));
  if(!f_mva) {
    printf("MVA file not found, exiting\n");
    return NULL;
  }

  TCanvas* c = new TCanvas("c","c",1000,600);
  TH1F *hMLP = (TH1F*) f_mva->Get("dpf");
  hMLP->SetName("hMLP");
  TH1F *hTot = (TH1F*) f_tot->Get("dpf");
  hTot->SetName("hTot");
  TH1F *hCuts = (TH1F*) hMLP->Clone("hCuts");
  hCuts->SetName("hCuts");

  TTree* t_cut = (TTree*) f_cut->Get("tmva_training_tree");
  if(!t_cut) {
    printf("Tree not found\n");
    return NULL;
  }
  printf("Drawing Cuts tree\n");

  t_cut->Draw("(p-pmc)>>hCuts", cuts);
  hCuts = (TH1F*) gDirectory->Get("hCuts");

  if(!hCuts || !hMLP) {
    printf("Histograms not found, exiting\n");
    return NULL;
  }

  gStyle->SetOptStat(1111);

  hCuts->SetLineWidth(2);
  hCuts->SetLineColor(kBlue);
  hCuts->SetFillColor(kAzure-4);
  hCuts->Draw("hist");

  hMLP->SetLineWidth(2);
  hMLP->SetLineColor(kRed);
  hMLP->SetFillStyle(3002);
  hMLP->SetFillColor(kRed);
  hMLP->Draw("hist sames");

  hTot->SetLineWidth(2);
  hTot->SetLineColor(kRed+3);
  // hTot->SetFillStyle(3002);
  // hTot->SetFillColor(kYellow);
  hTot->Draw("hist sames");

  hCuts->SetTitle("Track P - P (MC) Distribution");
  hCuts->SetXTitle("Reconstructed P - MC P (MeV)"); 
  hCuts->SetYTitle("Entries / .1 MeV");
  if(hCuts->GetMaximum() < hTot->GetMaximum())
    hCuts->SetMaximum(hTot->GetMaximum()*1.1);

  TF1* fCut, *fMLP;
  TLegend* legend = new TLegend();
  if(fitCrystal > 0) {
    printf("Declaring new fit functions\n");
    fCut = new TF1("fCut","crystal",-5.,5.,7); //fixing the (p-pmc) fitting range
    fCut->SetParameter(0,1); //norm
    fCut->SetParameter(1,0.); //mean
    fCut->SetParameter(2,1); //sigma
    fCut->SetParameter(3,1); //alpha Neg
    fCut->SetParameter(4,5); //power Neg
    fCut->SetParameter(5,1); //alpha Pos
    fCut->SetParameter(6,5); //power Pos
    fCut->SetParLimits(0,1,1e6);
    fCut->SetParLimits(1,-1.,1.);
    fCut->SetParLimits(2,0.01,1.);
    fCut->SetParLimits(3,0.1,10);
    fCut->SetParLimits(4,1,100);
    fCut->SetParLimits(5,0.1,10);
    fCut->SetParLimits(6,1,100);
    fMLP = new TF1("fMLP","crystal",-5.,5.,7); //fixing the (p-pmc) fitting range
    fMLP->SetParameter(0,1); //norm
    fMLP->SetParameter(1,0.); //mean
    fMLP->SetParameter(2,1); //sigma
    fMLP->SetParameter(3,1); //alpha Neg
    fMLP->SetParameter(4,5); //power Neg
    fMLP->SetParameter(5,1); //alpha Pos
    fMLP->SetParameter(6,5); //power Pos
    fMLP->SetParLimits(0,1,1e6);
    fMLP->SetParLimits(1,-1.,1.);
    fMLP->SetParLimits(2,0.01,1.);
    fMLP->SetParLimits(3,0.1,10);
    fMLP->SetParLimits(4,1,100);
    fMLP->SetParLimits(5,0.1,10);
    fMLP->SetParLimits(6,1,100);
    printf("Fitting functions\n");
    TFitResultPtr r1 = hCuts->Fit(fCut,"RS");
    TFitResultPtr r2 = hMLP->Fit(fMLP,"RS");
    fCut->Draw("same");
    fMLP->Draw("same");
    fCut->SetLineColor(kAzure+5);
    fCut->SetLineWidth(3);
    fMLP->SetLineColor(kOrange+1);
    fMLP->SetLineWidth(3);
    gStyle->SetOptFit(0);
    legend->AddEntry(hCuts,"Cuts");
    legend->AddEntry(hMLP,"MLP");
    legend->AddEntry(hTot,"Uncut");
    legend->AddEntry(fCut,Form("Cut Crystal #sigma = %.3f PPos = %.3f",r1->Parameter(2), r1->Parameter(6)));
    legend->AddEntry(fMLP,Form("MLP Crystal #sigma = %.3f PPos = %.3f",r2->Parameter(2), r2->Parameter(6)));

  }  else {
    legend->AddEntry(hCuts,"Cuts");
    legend->AddEntry(hMLP,"MLP");
    legend->AddEntry(hTot,"Uncut");
  }

  if(logPlot)  c->SetLogy();
  c->SetGridx();
  c->SetGridy();
  // c->SetFillColor(kGray);
  // c->GetFrame()->SetFillColor(kWhite);

  legend->Draw();
  printf("Efficiency MLP = %.3f\nEfficiency Cut = %.3f\n", hMLP->GetEntries()*1./hTot->GetEntries(),
	 hCuts->GetEntries()*1./hTot->GetEntries());
  c->Print(Form("dpf_%s_bkgweight_%i_plot.png",Signal,bkg_weight));
  return c;
}

int make_plots(const char* Algorithm = "DAR", const char* Signal = "ep") {

  TString tmvaName, alg, sig;

  alg = Algorithm;
  alg.ToLower();
  sig = Signal;
  sig.ToLower();
  
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
    return -1;
  }
  return 0;
}
