// [[Rcpp::depends(fntl)]]
#include <Rcpp.h>
#include "fntl.h"

#include "metaprediction.h"

using namespace Rcpp;
using namespace std;


// Edgington combined p-value function
// [[Rcpp::export]]
Rcpp::NumericVector pfct_edge_cpp(Rcpp::NumericVector h0,
                                  Rcpp::NumericVector es,
                                  Rcpp::NumericVector se) {
  
  
  // Number of studies
  int k = es.size();
  
  // Evaluated h0
  int Lh0 = h0.size();
  
  // Matrix to store individual study p-values
  NumericMatrix ps(Lh0, k);
  
  // For each study, compute one-sided Wald p-values
  for (int j = 0; j < k; ++j) {
    for (int i = 0; i < Lh0; ++i) {
      ps(i, j) = 1.0 - R::pnorm(es[j], h0[i], se[j], 1, 0);
    }
  }
  
  // Vector to store Edgington combined p-values
  Rcpp::NumericVector pcombined(Lh0);
  
  // Compute Edgington combined p-values
  for (int i = 0; i < Lh0; ++i) {
    Rcpp::NumericVector row(k);
    
    // P-value of each study for a specific H0
    for (int j = 0; j < k; ++j) {
      row[j] = ps(i, j);
    }
    
    // Edgington combination (with normal approx for k >= 12)
    double s = sum(row);
    double pE = 0.0;
    
    if (k < 12) {
      int sumUpper = floor(s);
      for (int ll = 0; ll <= sumUpper; ++ll) {
        pE += std::pow(-1, ll) * R::choose(k, ll) * std::pow(s - ll, k);
      }
      pE /= R::gammafn(k + 1);
    } else {
      double z = std::sqrt(12.0 * k) * (s / k - 0.5);
      pE = R::pnorm(z, 0.0, 1.0, 1, 0);
    }
    
    pcombined[i] = pE;
  }
  
  return pcombined;
}

// One-sided Wald p-value function
// [[Rcpp::export]]
Rcpp::NumericVector p_wald(double x,
                           Rcpp::NumericVector es,
                           Rcpp::NumericVector se) {
  
  Rcpp::NumericVector p = es.size();
  
  for (int ii = 0; ii < es.size(); ++ii) {
    p[ii] = R::pnorm(es[ii], x, se[ii], false, false);
  }
  
  return p;
}

// Optimize Edgington p-value function
// [[Rcpp::export]]
double opti_edge(Rcpp::NumericVector es,
                 Rcpp::NumericVector se) {
  
  // Find mean(p) = 0.5
  fntl::dfd f = [&](double x) {
    
    // One-sided p-values
    Rcpp::NumericVector p = p_wald(x, es, se);
    
    // Sum of p-values
    double s = 0.0;
    
    for (int ii = 0; ii < es.size(); ++ii) {
      s += p[ii];
    }
    
    // Mean of p-values
    s/=es.size();
    
    return s - 0.5;
    
  };
  
  // Root-finding settings
  fntl::findroot_args args;
  args.action = fntl::error_action::NONE;
  
  // Boundaries for root-finding
  double min_es = *std::min_element(es.begin(), es.end()) - *std::max_element(se.begin(), se.end());
  double max_es = *std::max_element(es.begin(), es.end()) + *std::max_element(se.begin(), se.end());
  
  // Find the root using Brent algorithm
  auto out = fntl::findroot_brent(f, min_es, max_es, args);
  
  // Return
  return out.root;
  
}

// Confidence density of mu (evaluated for one mu)
// [[Rcpp::export]]

double cd_single(Rcpp::NumericVector h0,
                 Rcpp::NumericVector es,
                 Rcpp::NumericVector se,
                 double h = 1e-4) {
  
  // Function for which derivative is computed
  fntl::dfv f = [&](Rcpp::NumericVector x) {
    return pfct_edge_cpp(x, es, se)[0];
  };
  
  // Compute and return the numerical derivative
  return fntl::fd_deriv(f, h0, 0, h);
}

// Confidence density of mu (evaluates multiple mu's)
// [[Rcpp::export]]
Rcpp::NumericVector CD_cpp(Rcpp::NumericVector h0,
                           Rcpp::NumericVector es,
                           Rcpp::NumericVector se,
                           double h = 1e-4) { // same as in numDeriv::grad
  
  Rcpp::NumericVector dv(h0.size());
  
  for (int ii = 0; ii < h0.size(); ++ii) {
    dv[ii] = cd_single(Rcpp::NumericVector::create(h0[ii]), es, se, h);
  }
  
  return dv;
  
}











