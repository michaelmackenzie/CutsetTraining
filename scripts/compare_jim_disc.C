//Script to compare Feldman-Cousins discovery to Jim's
Significances sig_;

void compare_jim_disc() {
  double data[] = {
    //nbkg   nsig       SES           upper limit      N(discover)      R discovery
    1.336,   4.225,   2.37E-17,        8.32E-17,        1.10E+01,        2.60E-16,
    1.186,   4.15 ,   2.41E-17,        8.22E-17,        1.00E+01,        2.41E-16,
    1.046,   4.071,   2.46E-17,        8.13E-17,        1.00E+01,        2.46E-16,
    0.933,   3.983,   2.51E-17,        8.09E-17,        1.00E+01,        2.51E-16,
    0.842,   3.883,   2.58E-17,        8.11E-17,        9.00E+00,        2.32E-16,
    0.768,   3.779,   2.65E-17,        8.17E-17,        9.00E+00,        2.38E-16,
    0.7  ,   3.67 ,   2.72E-17,        8.25E-17,        9.00E+00,        2.45E-16,
    0.641,   3.546,   2.82E-17,        8.38E-17,        8.00E+00,        2.26E-16,
    0.583,   3.415,   2.93E-17,        8.55E-17,        8.00E+00,        2.34E-16,
    0.544,   3.27 ,   3.06E-17,        8.84E-17,        8.00E+00,        2.45E-16,
    0.508,   3.12 ,   3.20E-17,        9.17E-17,        8.00E+00,        2.56E-16,
    0.475,   2.957,   3.38E-17,        9.59E-17,        8.00E+00,        2.71E-16,
    0.451,   2.789,   3.59E-17,        1.01E-16,        7.00E+00,        2.51E-16,
    0.43 ,   2.608,   3.83E-17,        1.07E-16,        7.00E+00,        2.68E-16,
    0.409,   2.422,   4.13E-17,        1.15E-16,        7.00E+00,        2.89E-16,
    0.374,   2.23 ,   4.48E-17,        1.24E-16,        7.00E+00,        3.14E-16,
    0.356,   2.029,   4.93E-17,        1.35E-16,        7.00E+00,        3.45E-16,
    0.341,   1.819,   5.50E-17,        1.50E-16,        7.00E+00,        3.85E-16,
    0.322,   1.608,   6.22E-17,        1.69E-16,        7.00E+00,        4.35E-16,
    0.308,   1.403,   7.13E-17,        1.93E-16,        7.00E+00,        4.99E-16,
    0.292,   1.203,   8.31E-17,        2.24E-16,        7.00E+00,        5.82E-16,
    0.278,   1.008,   9.92E-17,        2.66E-16,        6.00E+00,        5.95E-16,
    0.263,   0.821,   1.22E-16,        3.25E-16,        6.00E+00,        7.30E-16,
    0.249,   0.657,   1.52E-16,        4.04E-16,        6.00E+00,        9.13E-16,
    0.235,   0.515,   1.94E-16,        5.13E-16,        6.00E+00,        1.17E-15,
    0.222,   0.391,   2.56E-16,        6.73E-16,        6.00E+00,        1.53E-15,
    0.208,   0.293,   3.41E-16,        8.94E-16,        6.00E+00,        2.05E-15,
    0.195,   0.217,   4.60E-16,        1.20E-15,        6.00E+00,        2.76E-15,
    0.182,   0.161,   6.19E-16,        1.61E-15,        6.00E+00,        3.72E-15,
    0.169,   0.119,   8.37E-16,        2.17E-15,        6.00E+00,        5.02E-15,
    0.156,   0.09 ,   1.11E-15,        2.87E-15,        6.00E+00,        6.69E-15,
    0.143,   0.067,   1.48E-15,        3.81E-15,        5.00E+00,        7.42E-15,
    0.13 ,   0.051,   1.95E-15,        4.99E-15,        5.00E+00,        9.76E-15,
    0.117,   0.038,   2.62E-15,        6.67E-15,        5.00E+00,        1.31E-14,
    0.104,   0.028,   3.52E-15,        8.90E-15,        5.00E+00,        1.76E-14,
    0.092,   0.021,   4.66E-15,        1.17E-14,        5.00E+00,        2.33E-14,
    0.079,   0.016,   6.26E-15,        1.57E-14,        5.00E+00,        3.13E-14,
    0.066,   0.011,   8.90E-15,        2.22E-14,        5.00E+00,        4.45E-14,
    0.053,   0.008,   1.22E-14,        3.04E-14,        4.00E+00,        4.89E-14,
    0.04 ,   0.006,   1.71E-14,        4.22E-14,        4.00E+00,        6.83E-14,
    0.026,   0.004,   2.70E-14,        6.64E-14,        4.00E+00,        1.08E-13,
    0.013,   0.002,   5.22E-14,        1.28E-13,        3.00E+00,        1.57E-13
  };
  int ndata = sizeof(data)/sizeof(*data)/6;
  TH1F* hdiffmedUL = new TH1F("hdiffmedul", "Relative median UL difference", 40., -0.3, 0.3);
  TH1F* hdiffavgUL = new TH1F("hdiffavgul", "Relative mean UL difference", 40., -0.3, 0.3);
  TH1F* hdiffmed5S = new TH1F("hdiffmed5s", "Relative median discovery difference", 40., -0.3, 0.3);
  TH1F* hdiffavg5S = new TH1F("hdiffavg5s", "Relative mean discovery difference", 40., -0.3, 0.3);
  for(int idata = 0; idata < ndata; ++idata) {
    double nbkg  = data[6*idata  ];
    double nsig  = data[6*idata+1];
    double ses   = data[6*idata+2];
    double ul    = data[6*idata+3];
    int    ndisc = data[6*idata+4];
    double rdisc = data[6*idata+5];
    int NDISC = sig_.FCSigmaN(nbkg, 5.);
    if(ndisc != NDISC) {
      printf("N(discovery) disagrees! J = %i, M = %i, nbkg = %.3e\n", ndisc, NDISC, nbkg);
    }
    double UL_MED = ses*sig_.FCMedianUpperLimit(nbkg, 0.9);
    double UL_AVG = ses*sig_.FCUpperLimit(nbkg, 0.9);
    double DC_MED = ses/sig_.FCDiscovery(1., nbkg);
    double DC_AVG = ses/sig_.AverageFiveSigma(1., nbkg);
    hdiffmedUL->Fill((ul - UL_MED)/(ul+UL_MED)/2.);
    hdiffavgUL->Fill((ul - UL_AVG)/(ul+UL_AVG)/2.);
    hdiffmed5S->Fill((rdisc - DC_MED)/(rdisc+DC_MED)/2.);
    hdiffavg5S->Fill((rdisc - DC_AVG)/(rdisc+DC_AVG)/2.);
  }
  TCanvas* c = new TCanvas("c_compare", "c_compare", 1200, 1000);
  c->Divide(2,2);
  c->cd(1);
  hdiffmedUL->Draw();
  hdiffmedUL->SetLineWidth(2);
  hdiffmedUL->SetXTitle("#Delta UL / Avg UL");
  c->cd(2);
  hdiffavgUL->Draw();
  hdiffavgUL->SetLineWidth(2);
  hdiffavgUL->SetXTitle("#Delta UL / Avg UL");
  c->cd(3);
  hdiffmed5S->Draw();
  hdiffmed5S->SetLineWidth(2);
  hdiffmed5S->SetXTitle("#Delta R(discovery) / Avg R(discovery)");
  c->cd(4);
  hdiffavg5S->Draw();
  hdiffavg5S->SetLineWidth(2);
  hdiffavg5S->SetXTitle("#Delta R(discovery) / Avg R(discovery)");
  c->Print("figures/compare_jim_disc.png");
}
