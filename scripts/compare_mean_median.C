//compare FC discovery method
TCanvas* compare_mean_median(double nsignal, double nbackground, bool print = false) {

  gStyle->SetOptStat(0);
  int nmax = nsignal + nbackground + 5.*sqrt(nbackground) + 5.*sqrt(nsignal+nbackground);
  TH1F* hbkg = new TH1F("hbkg", "hbkg", nmax+1, 0, nmax+1);
  TH1F* hsig = new TH1F("hsig", "hsig", nmax+1, 0, nmax+1);
  TH1F* hbkgcdf = new TH1F("hbkgcdf", "hbkgcdf", nmax+1, 0, nmax+1);
  TH1F* hbkgsig = new TH1F("hbkgsig", "hbkgsig", nmax+1, 0, nmax+1);
  TH1F* hsigma = new TH1F("hsigma", "hsigma", 100, -5., 10.);
  double median(-1.), prob(0.), mean(0.);
  for(int nseen = 0; nseen <= nmax; ++nseen) {
    const double p_bkg = ROOT::Math::poisson_pdf(nseen, nbackground);
    const double p_sig = ROOT::Math::poisson_pdf(nseen, nbackground+nsignal);
    //if N seen below background mean, flip cdf direction
    const double p_bkgsig = (nseen < nbackground) ? ROOT::Math::poisson_cdf(nseen, nbackground) : 1. - ROOT::Math::poisson_cdf(nseen-1, nbackground);
    double sigma = ROOT::Math::gaussian_quantile(1. - p_bkgsig, 1.);
    if(isinf(sigma)) {
      sigma = (p_bkg > 0.5) ? -10. : 10.;
      cout << "!!! Sigma value is inf!\n";
    } else if(nseen < p_bkg) {
      sigma *= 1;
    }
    cout << "N(seen) = " << nseen << ", sigma = " << sigma << endl;
    prob += p_sig;
    if(prob >= 0.5 && median < 0.) median = sigma;
    mean += sigma*p_sig;
    hbkg->SetBinContent(nseen+1, p_bkg);
    hsig->SetBinContent(nseen+1, p_sig);
    hbkgcdf->SetBinContent(nseen+1, p_bkgsig);
    hbkgsig->SetBinContent(nseen+1, sigma);
    hsigma->Fill(sigma, p_sig);
  }

  TCanvas* c = new TCanvas("c_compare", "c_compare", 1400, 800);
  c->Divide(2,1);
  auto pad = c->cd(1);

  /////////////////////////////
  // Draw N(bkg) and N(sig)  //
  /////////////////////////////

  hsig->SetLineColor(kRed);
  hsig->SetLineWidth(2);
  hsig->SetFillStyle(3005);
  hsig->SetFillColor(hsig->GetLineColor());

  hbkg->SetLineWidth(2);
  hbkg->SetFillStyle(3004);
  hbkg->SetFillColor(hbkg->GetLineColor());
  // hbkgcdf->SetFillStyle(3004);
  // hbkgcdf->SetFillColor(hbkgcdf->GetLineColor());
  hbkgsig->SetFillStyle(3005);
  hbkgsig->SetFillColor(hbkgsig->GetLineColor());

  hbkg->Draw("hist");
  hsig->Draw("hist same");
  // hbkgcdf->Draw("hist same");
  // hbkgsig->Scale(1./hbkgsig->GetMaximum());
  // hbkgsig->Draw("hist same");
  hbkg->GetYaxis()->SetRangeUser(1.e-9, 10.);
  hbkg->SetTitle("Background and Signal PDFs");
  hbkg->SetXTitle("N(events)");
  hbkg->SetYTitle("P(N)");
  hbkg->GetYaxis()->SetTitleOffset(1.2);

  TLegend* leg = new TLegend(0.38, 0.76, 0.9, 0.9);
  leg->AddEntry(hbkg, Form("P(n | #mu_{b}),  #mu_{b} = %.2f", nbackground));
  leg->AddEntry(hsig, Form("P(n | #mu_{b}+#mu_{s}), #mu_{s} = %.2f", nsignal));
  leg->SetTextSize(0.04);
  leg->Draw();
  pad->SetLogy();

  /////////////////////////////
  // Draw sigma distribution //
  /////////////////////////////

  c->cd(2);
  hsigma->SetLineWidth(2);
  hsigma->SetFillColor(hsigma->GetLineColor());
  hsigma->SetFillStyle(3003);
  hsigma->Draw("hist");
  TLine* lmean = new TLine(mean, 0., mean, 1.05*hsigma->GetMaximum());
  lmean->SetLineWidth(2);
  lmean->SetLineColor(kRed);
  lmean->Draw("same");
  TLine* lmedian = new TLine(median, 0., median, 1.05*hsigma->GetMaximum());
  lmedian->SetLineWidth(2);
  lmedian->SetLineColor(kGreen);
  lmedian->Draw("same");

  hsigma->SetTitle("#sigma deviation from H(background only)");
  hsigma->SetXTitle("#sigma deviation");
  hsigma->SetYTitle("P(#sigma)");
  hsigma->GetYaxis()->SetTitleOffset(1.5);
  hsigma->GetYaxis()->SetRangeUser(0., 1.3*hsigma->GetMaximum());

  TLegend* leg2 = new TLegend(0.1, 0.75, 0.6, 0.9);
  leg2->AddEntry(hsigma, "#sigma deviation");
  leg2->AddEntry(lmean, Form("Mean #sigma (#bar{#sigma}) = %.3f", mean), "L");
  leg2->AddEntry(lmedian, Form("Median #sigma = %.3f", median), "L");
  leg2->Draw();

  if(print)
    c->Print("figures/compare_mean_median.png");
  return c;
}
