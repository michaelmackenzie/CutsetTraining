// Compare FC discovery and 90% CL for signal and background numbers

Significances* sigs_;
int verbose_ = 0;
bool doAllModes_ = true;

void FC_Discovery_vs_Limit(double nsignal, double n_min, double n_max,
                               int nsteps = 1000, bool print = false, int Mode = 0) {

  sigs_ = new Significances();
  // sigs_->verbose_ = 2;
  double nbkgs[nsteps+1];
  double mudiscover[nsteps+1]; //discovery signal strength
  double mulimit[nsteps+1]; //average upper limit signal strength
  double medianlimit[nsteps+1]; //median upper limit signal strength
  double avgsig[nsteps+1]; //average discovery for example signal strength
  double muavg[nsteps+1]; //mu so average significance is 5 sigma
  double ndiscover[nsteps+1];
  double nlim[nsteps+1]; //median n seen for 90% CL
  double disc_scale = 1.;//(Mode < 0) ? 1. : 0.4;
  for(int istep = 0; istep <= nsteps; ++istep) {
    double nbkg = n_min + istep*(n_max - n_min)/nsteps;
    mudiscover[istep] = disc_scale*nsignal*1./sigs_->FCDiscovery(nsignal, nbkg);
    mulimit[istep] = sigs_->FCUpperLimit(nbkg, 0.9);
    int ilim;
    medianlimit[istep] = sigs_->FCMedianUpperLimit(nbkg, 0.9, ilim);
    nlim[istep] = ilim;
    ndiscover[istep] = sigs_->FCSigmaN(nbkg, 5.);
    avgsig[istep] = sigs_->AverageSignificance(nsignal,nbkg);
    muavg[istep] = disc_scale*nsignal/sigs_->AverageFiveSigma(nsignal,nbkg);
    nbkgs[istep] = nbkg;
    if(verbose_ > 0)
      printf("%i: bkg = %.3e, discover = %.3e, limit = %.3e, med. limit = %.3e, nmedlim = %i, ndiscover = %.0f, avgsig = %.3e\n",
             istep, nbkg, mudiscover[istep]/disc_scale, mulimit[istep], medianlimit[istep], ilim, ndiscover[istep], avgsig[istep]);
  }
  TGraph* gdiscover = new TGraph(nsteps+1, nbkgs, mudiscover);
  gdiscover->SetName("gdiscover");
  gdiscover->SetLineColor(kBlue);
  gdiscover->SetMarkerColor(kBlue);
  gdiscover->SetLineWidth(2);
  gdiscover->SetMarkerStyle(20);
  gdiscover->SetMarkerSize(0.8);

  TGraph* glimit = new TGraph(nsteps+1, nbkgs, mulimit);
  glimit->SetName("glimit");
  glimit->SetLineColor(kRed);
  glimit->SetMarkerColor(kRed);
  glimit->SetLineWidth(2);
  glimit->SetMarkerStyle(20);
  glimit->SetMarkerSize(0.8);

  TGraph* gmedlim = new TGraph(nsteps+1, nbkgs, medianlimit);
  gmedlim->SetName("gmedlim");
  gmedlim->SetLineColor(kOrange-3);
  gmedlim->SetMarkerColor(kOrange-3);
  gmedlim->SetLineWidth(2);
  gmedlim->SetMarkerStyle(20);
  gmedlim->SetMarkerSize(0.8);

  TGraph* gavgsig = new TGraph(nsteps+1, nbkgs, avgsig);
  gavgsig->SetName("gavgsig");
  gavgsig->SetLineColor(kOrange-3);
  gavgsig->SetMarkerColor(kOrange-3);
  gavgsig->SetLineWidth(2);
  gavgsig->SetMarkerStyle(20);
  gavgsig->SetMarkerSize(0.8);

  TGraph* gavgmu = new TGraph(nsteps+1, nbkgs, muavg);
  gavgmu->SetName("gavgmu");
  gavgmu->SetLineColor(kMagenta);
  gavgmu->SetMarkerColor(kMagenta);
  gavgmu->SetLineWidth(2);
  gavgmu->SetMarkerStyle(20);
  gavgmu->SetMarkerSize(0.8);

  if(doAllModes_) Mode = -1;
  do {
    if(Mode >= 0)
      gdiscover->SetTitle("Discovery and Upper Limit vs Background;Background #mu;Signal #mu");
    else
      gdiscover->SetTitle("Discovery vs Background;Background #mu;Signal #mu");

    TCanvas* c = new TCanvas(Form("c_fc_%i", Mode), Form("c_fc_%i", Mode), 1000, 700);
    if(Mode != 4) {
      gdiscover->Draw("AP");
      gdiscover->GetYaxis()->SetRangeUser(0.01, mudiscover[nsteps]*1.1);
      gdiscover->GetXaxis()->SetRangeUser(n_min, n_max);
      c->Update();
    } else {
      glimit->Draw("AP");
      glimit->GetYaxis()->SetRangeUser(0., mudiscover[nsteps]*1.1);
      glimit->GetXaxis()->SetRangeUser(n_min, n_max);
      glimit->SetTitle("Upper Limit vs Background;Background #mu;Signal #mu");
      c->Update();
    }

    if(Mode == 1) {
      gavgsig->Draw("P");
    }
    if(Mode == 1 || Mode == 3) {
      gavgmu->Draw("P");
    }
    if(Mode == 0 || Mode == 3) {
      glimit->Draw("P");
    }
    if(Mode == 2 || Mode == 3 || Mode == 4) {
      gmedlim->Draw("P");
    }

    if(Mode >= 0) {
      TLegend* leg = new TLegend(0.1, 0.7, 0.35, 0.9);
      if(Mode != 4)
        leg->AddEntry(gdiscover, Form("#mu(#tilde{#sigma} #geq 5)%s", (disc_scale < 1.) ? Form("*%.1f", disc_scale) : "" ), "L");
      if(Mode == 1 || Mode == 3)
        leg->AddEntry(gavgmu, Form("#mu(#bar{#sigma} = 5)%s", (disc_scale < 1.) ? Form("*%.1f", disc_scale) : "" ), "L");
      if(Mode == 2 || Mode == 3 || Mode == 4)
        leg->AddEntry(gmedlim, "#tilde{#mu}(90% CL)", "L");
      if(Mode == 0 || Mode == 3 || Mode == 4)
        leg->AddEntry(glimit, "#bar{#mu}(90% CL)", "L");
      if(Mode == 1)
        leg->AddEntry(gavgsig, Form("#bar{#sigma}(#mu=%.1f)", nsignal), "L");
      leg->Draw();
    }

    //add a second axis for the N(Five Sigma)
    double rightmax = (Mode == 4) ? 1.1*nlim[nsteps] : 1.1*ndiscover[nsteps];
    double naxis[nsteps+1];
    for(int istep = 0; istep <= nsteps; ++istep) {
      if(Mode == 4)
        naxis[istep] = nlim[istep]*gPad->GetUymax()/rightmax;
      else
        naxis[istep] = ndiscover[istep]*gPad->GetUymax()/rightmax;
    }

    TGaxis* axis = new TGaxis(gPad->GetUxmax(), gPad->GetUymin(),
                              gPad->GetUxmax(), gPad->GetUymax(), 0, rightmax,
                              510, "+L");
    axis->SetLineColor (kGreen+2);
    axis->SetLabelColor(kGreen+2);
    axis->SetTitleColor(kGreen+2);
    axis->SetTitle((Mode == 4) ? "Median N(CL)" : "N(5#sigma)");
    axis->Draw();

    TGraph* gndisc = new TGraph(nsteps+1, nbkgs, naxis);
    gndisc->SetName(Form("gndisc_%i", Mode));
    gndisc->SetLineColor(kGreen+2);
    gndisc->SetLineWidth(2);
    gndisc->SetMarkerColor(kGreen+2);
    gndisc->SetMarkerStyle(20);
    gndisc->SetMarkerSize(0.8);
    gndisc->Draw("P");

    if(print) {
      c->Print(Form("figures/FC_Discovery_vs_Limit_Mode_%i.png", Mode));
    }
    ++Mode;
  } while(Mode < 5 && doAllModes_);

}
