#ifndef __Significances__HH
#define __Significances__HH

#include "Math/ProbFuncMathCore.h"
#include "Math/PdfFuncMathCore.h"
#include "Math/QuantFuncMathCore.h"
#include "TRandom3.h"
#include <set>
class Significances : public TObject {
public:


public:

  Significances(int seed = 90): maxAttempts_(1e3), seed_(seed), verbose_(0) {
    rnd_ = new TRandom3(seed_);
  }

  //Feldman-Cousins confidence belt for a signal mu, from n_min to n_max at confidence level cl
  void FCBeltMu(double nsignal, double nbackground, double cl, int& n_min, int& n_max) {
    double prob = 0.;
    n_min = 1e9; n_max = 0;
    if(cl >= 1. || cl <= 0. || nsignal < 0. || nbackground < 0. || (nsignal == 0. && nbackground == 0.)) return;
    std::set<int> nfound; //set of ns added to the probability, use a set so search is O(1)
    const double nsigma = ROOT::Math::gaussian_quantile(cl, 1);
    const double mutot = nsignal + nbackground;
    unsigned nmax = 5. + (mutot) + (nsigma+3.)*sqrt(mutot);
    while(prob < cl) { //continue adding n's to the band until the probibility matches the confidence level
      //find the next n with the highest F-C ordering parameter, R = P(n | mu) / P(n | max(B, n - B))
      if(nfound.size() <= nmax) nmax += 10; //ensure there are always ns to check
      double rmax = -1.;
      int currn = -1;
      for(unsigned n = 0; n < nmax; ++n) {
        if(nfound.find(n) != nfound.end()) continue;
        double r = ROOT::Math::poisson_pdf(n, mutot) / ROOT::Math::poisson_pdf(n, std::max((double) n, nbackground));
        if(r > rmax) {rmax = r; currn = n;}
      }
      nfound.insert(currn);
      prob += ROOT::Math::poisson_pdf(currn, mutot);
      n_min = std::min(n_min, currn);
      n_max = std::max(n_max, currn);
    }
  }

  //Feldman-Cousins based median upper limit
  double FCMedianUpperLimit(double mu_background, double cl, int& ncl) {
    ncl = -1;
    double prob = 0.;
    //find the value of n that crosses the 50% threshold --> this is the median
    while(prob < 0.5) {
      ++ncl;
      prob += ROOT::Math::poisson_pdf(ncl, mu_background);
    }
    double limit = GetFCUpperLimit(ncl, mu_background, cl);
    return limit;
  }

  //Feldman-Cousins based median upper limit
  double FCMedianUpperLimit(double mu_background, double cl = 0.9) {
    int ncl;
    return FCMedianUpperLimit(mu_background, cl, ncl);
  }

  //Average Feldman-Cousins based upper limit
  double FCUpperLimit(double mu_background, double cl = 0.9) {
    const static int nsteps = 10 + mu_background + (ROOT::Math::gaussian_quantile(cl, 1) + 2.)*sqrt(mu_background);
    double limit_val = 0.;
    //get the upper limit for a range of N's, weighted by their probability
    for(int istep = 0; istep < nsteps; ++istep) {
      double limit_step = GetFCUpperLimit(istep, mu_background, cl);
      double prob = ROOT::Math::poisson_pdf(istep, mu_background);
      limit_val += limit_step * prob;
    }
    return limit_val;
  }

  // //mu if background is high
  // double HighNUpperLimit(double mu_background, double cl = 0.9) {
  //   const static int nattempts = 500;
  //   double limit_val = 0.;
  //   //get the upper limit for a range of N's, weighted by their probability
  //   for(int attempt = 0; attempt < nattempts; ++attempt) {
  //     int n = rnd_->Poisson(mu_background);
  //     double limit_step = GetUpperLimitMu(n, cl) - mu_background;
  //     limit_val += limit_step;
  //   }
  //   limit_val /= nattempts;
  //   return limit_val;
  // }

  //for a given observation and background expectation, get FC upper limit
  double GetFCUpperLimit(int nseen, double mu_background, double cl) {
    //find the largest value of signal mu that this n is contained within its confidence belt
    const static double tolerance = 0.001; //mu_signal tolerance
    double mu_max = nseen + 10 + (ROOT::Math::gaussian_quantile(cl, 1) + 2.)*sqrt(nseen);
    double mu_min = 0.;
    double mu_curr = (mu_max + mu_min)/2.;
    while(abs(mu_max - mu_min) > tolerance) { //find the value of mu up to a tolerance
      int n_min, n_max;
      FCBeltMu(mu_curr, mu_background, cl, n_min, n_max);
      if(n_min <= nseen) mu_min = mu_curr; //mu includes median, so increase mu
      else             mu_max = mu_curr; //mu doesn't include median, so decrease
      mu_curr = (mu_max + mu_min)/2.;
    }
    return mu_curr;
  }

  //upper limit on n_signal if saw n_seen with an expectation of mu, not with Feldman-Cousins
  double GetUpperLimitMu(int n_seen, double cl = 0.9) {
    double limit = 1. - cl; //probability goal
    //starting guess for mu and the probability of it
    double mu_guess = std::max(n_seen + 1.7*std::sqrt(n_seen), 1.);
    double val = ROOT::Math::poisson_cdf(n_seen, mu_guess);

    //how close a value can be to the desired probability
    const static double tolerance = 1.e-3;
    const static double alpha = 0.1;
    int attempts = 0;

    while(abs(val-limit) > tolerance && attempts < maxAttempts_) {
      ++attempts;
      double change = mu_guess*alpha*(val-limit)/limit;
      mu_guess += change;
      val = ROOT::Math::poisson_cdf(n_seen, mu_guess);
    }
    if(attempts == maxAttempts_ && verbose_ > 0) {
      printf("!!! %s: Reached maximum attempts getting the upper limit value! mu_guess = %.3e, cl = %.3e, int n_seen = %i, val = %.3e\n",
             __func__, mu_guess, cl, n_seen, val);
    }
    return mu_guess;
  }

  //Number needed to see to be above the Feldman-Cousins one-sided n-sigma p-value confidence belt
  int FCSigmaN(const double nbackground, const double nsigma = 5.) {
    if(nbackground <= 0.) return 1; //need to see at least 1 event no matter what
    const double p_sigma = ROOT::Math::gaussian_cdf(nsigma, 1., 0.); //one-sided test
    int n_min, n_max;
    FCBeltMu(0., nbackground, p_sigma, n_min, n_max); //get the n sigma belt for the background only hypothesis
    return n_max+1; //number need to see to be beyond the n sigma belt for the background
  }

  //Feldman-Cousins based discovery
  double FCDiscovery(double nsignal, double nbackground) {
    if(nsignal < 0. || nbackground < 0.) {
      printf("!!! %s: Negative background or signal, returning -1!\n", __func__);
      return -1.;
    }
    double scale = 1.;
    const static double alpha = 0.1;
    const int nfive = FCSigmaN(nbackground, 5.); //number needed to see for five sigma discovery

    //find signal scale such that 50% claim a discovery
    int attempts = 0;
    double tolerance = 0.001;
    double val1 = 1. - ROOT::Math::poisson_cdf(nfive-1, nsignal*scale+nbackground); //probability of seeing >= nseen
    while(abs(val1 - 0.5) > tolerance && attempts < maxAttempts_) {
      ++attempts;
      scale *= 1. + (0.5-val1)*alpha;
      val1 = 1. - ROOT::Math::poisson_cdf(nfive-1, nsignal*scale+nbackground);
      if(verbose_ > 2) printf("--- %s: Loop %i for nfive, val1 = %.3e, scale = %.3e\n",
                              __func__, attempts, val1, scale);
      if(verbose_ > 0 && attempts == maxAttempts_)
        printf("!!! %s: Reached maximum attempts trying to find 5 sigma value! val1= %.3e, scale = %.3e\n",
               __func__, val1, scale);
    }
    return 1./scale;
  }

  //find the weighted average significance seen
  double AverageSignificance(double nsignal, double nbackground, const bool allow_negative = true) {
    if(nsignal < 0. || nbackground < 0.) {
      printf("!!! %s: Negative background or signal, returning -1!\n", __func__);
      return -1.;
    }
    const int nmax = nsignal + nbackground + 5.*sqrt(nbackground) + 3.*sqrt(nsignal);
    double avg = 0.;
    const static double pmin = 1.e-6;
    for(int nseen = 0; nseen <= nmax; ++nseen) {
      //probability see this or more in background only
      const double bkg_p = (nseen < nbackground) ? ROOT::Math::poisson_cdf(nseen, nbackground) : 1. - ROOT::Math::poisson_cdf(nseen - 1, nbackground);
      //get significance of this nseen
      double sigma = ROOT::Math::gaussian_quantile(1. - bkg_p, 1.);
      if(!allow_negative) sigma = std::max(0., sigma);
      if(isinf(sigma)) {
        if(verbose_ > 1) {
          printf("!!! %s: Sigma value is inf for mu_bkg = %.3e and nseen = %i\n", __func__, nbackground, nseen);
        }
        sigma = (bkg_p > 0.5) ? -10. : 10.; //figure out if extremely likely or unlikely, use 10 sigma as max
      }
      //probability of this for signal + background
      const double sig_p = ROOT::Math::poisson_pdf(nseen, nbackground+nsignal);
      //add this significance weighted by its probability
      avg += sigma*sig_p;
      if(verbose_ > 2) printf("--- %s: N(seen) = %i P(bkg) = %.3e --> Prob(%.3f sigma) = %.3e\n", __func__, nseen, bkg_p, sigma, sig_p);
      if(bkg_p > 1. - pmin && sig_p < pmin) break;
    }
    return avg;
  }

  //find signal strength that has average significance of 5 sigma
  double AverageFiveSigma(double nsignal, double nbackground, const bool allow_negative = true) {
    if(nsignal < 0. || nbackground < 0.) {
      printf("!!! %s: Negative background or signal, returning -1!\n", __func__);
      return -1.;
    }
    double scale = (nsignal+nbackground)/FCSigmaN(nbackground, 5.);
    const double alpha = 0.1;
    const double tolerance = 1.e-3;
    int attempts = 0;
    while(attempts < maxAttempts_) {
      ++attempts;
      double nsigma = AverageSignificance(scale*nsignal, nbackground, allow_negative);
      if(verbose_ > 1) {
        printf("--- %s: Attempt %i has average significance %.3f for nsignal*scale = %.3e\n",
               __func__, attempts, nsigma, scale*nsignal);
      }
      if(abs(nsigma - 5.) < tolerance) break;
      if(attempts == maxAttempts_) {
        printf("!!! %s: Hit maximum attempts, returning sigma = %.3e scaling!\n", __func__, nsigma);
      } else {
        scale *= std::min(2.,std::max(0.2, 1. + alpha*(5.-nsigma)/(5.)));
      }
    }
    return 1./scale; //optimizers maximize, not minimize
  }

  int maxAttempts_;
  int seed_;
  int verbose_;

  TRandom3* rnd_;
};

#endif
