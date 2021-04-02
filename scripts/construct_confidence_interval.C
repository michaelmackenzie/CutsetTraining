int verbose_ = 1;

void construct_confidence_interval(double mu_s, double mu_b, double alpha, double sigma, int& nmin_f, int& nmax_f) {
  if(sigma > 0.) {
    double p_onesided = ROOT::Math::gaussian_cdf(sigma, 1., 0.);
    alpha = 2.*p_onesided - 1.;
  }
  int n0 = mu_b+mu_s;
  int nmin = max( 0, (int) (n0 - 10.*sqrt(n0)));
  int nmax = max(10, (int) (n0 + 10.*sqrt(n0)));

  double prob = 0.;
  map<int, bool> ns;
  nmax_f = -1; nmin_f = nmax;
  while(prob < alpha) {
    double ratio(-1.);
    int nbest = -1;
    double pbest(0.);
    if(nmin_f == nmin && nmin > 0) nmin = max(nmin-10, 0);
    if(nmax_f == nmax) nmax += 10;
    for(int n = nmin; n <= nmax; ++n) {
      if(ns[n]) continue;
      double mu_best = max(mu_b, n - mu_b);
      double p = ROOT::Math::poisson_pdf(n, mu_s+mu_b);
      double pb = ROOT::Math::poisson_pdf(n, mu_best);
      if(ratio < p/pb) {ratio = p/pb; nbest = n; pbest=p;}
    }
    ns[nbest] = true;
    prob += pbest;
    if(verbose_ > 1) cout << "nbest = " << nbest << ", ratio = " << ratio << ", prob = " << prob << endl;
    if(nmax_f < nbest) nmax_f = nbest;
    if(nmin_f > nbest) nmin_f = nbest;
  }
  if(verbose_ > 0) {
    cout << "Final range = " << nmin_f << " - " << nmax_f << endl;
    cout << "Overcoverage = " << prob - alpha << endl;
  }
}

void construct_confidence_interval(double mu_s, double mu_b, double alpha = 0.9, double sigma = -1) {
  int nminf, nmaxf;
  construct_confidence_interval(mu_s, mu_b, alpha, sigma, nminf, nmaxf);
}

TCanvas* construct_confidence_interval() {
  double mu_s_min =  0.;
  double mu_s_max = 20.;
  int nsteps = 100;
  double mu_b = 0.4;
  double alpha = 0.9;
  double sigma = -1;
  double mus[nsteps+1], centers[nsteps+1], widths[nsteps+1], heights[nsteps+1];
  verbose_ = 0;
  int n_obs = 4;
  double b_obs_min(10*n_obs), b_obs_max(-1);
  for(int istep = 0; istep <= nsteps; ++istep) {
    int n1, n2;
    double mu_s = mu_s_min + istep*(mu_s_max - mu_s_min)/nsteps;
    construct_confidence_interval(mu_s, mu_b, alpha, sigma, n1, n2);
    if(n1 <= n_obs && n_obs <= n2) { //included in observed
      b_obs_min = min(b_obs_min, mu_s);
      b_obs_max = max(b_obs_max, mu_s);
    }
    centers[istep] = (n1+n2)/2.;
    widths[istep] = (n2-n1)/2.;
    heights[istep] = (mu_s_max-mu_s_min)/nsteps;
    mus[istep] = mu_s;
  }
  TGraph* g = new TGraphErrors(nsteps+1, centers, mus, widths, heights);
  g->SetName("g_belt");
  g->SetTitle(Form("Confidence belts for #alpha = %.2f and #mu_{b} = %.2f;n;#mu_{s}", alpha, mu_b));
  g->SetLineWidth(2);
  g->SetLineColor(kBlue);
  g->SetFillColor(kBlue);
  g->SetFillStyle(3001);
  TCanvas* c = new TCanvas("c_belt", "c_belt", 1000, 800);
  g->Draw("AE2");
  g->GetYaxis()->SetRangeUser(mu_s_min, mu_s_max);

  double nobs[] = {(double) n_obs+0.5}; double obsxe[] = {0.};
  double obsyc[] = {(b_obs_min+b_obs_max)/2.};
  double obsye[] = {(b_obs_max-b_obs_min)/2.};
  TGraph* gobs = new TGraphErrors(1, nobs, obsyc, obsxe, obsye);
  gobs->SetName("g_obs");
  gobs->SetLineColor(kRed);
  gobs->SetLineWidth(2);
  // gobs->Draw("E");

  c->Print("figures/confidence_belt.png");

  cout << "For an example " << n_obs << " events, the belt is " << b_obs_min << " - " << b_obs_max << endl;
  return c;
}
