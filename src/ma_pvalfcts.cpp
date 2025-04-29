// [[Rcpp::depends(fntl)]]
#include <Rcpp.h>
#include "fntl.h"
using namespace Rcpp;


// Edgingtons combined p-value function
// [[Rcpp::export]]
NumericVector pfct_edge_cpp(NumericVector h0,
                            NumericVector es,
                            NumericVector se) {


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
  NumericVector pcombined(Lh0);

  // Compute Edgington combined p-values
  for (int i = 0; i < Lh0; ++i) {
    NumericVector row(k);

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


// Function to generate on sample of mu from the confidence distribution
// [[Rcpp::export]]

double sample_one_mu(double s_tau2,
                     NumericVector es,
                     NumericVector se) {

  // Adjust the standard errors based on tau2
  NumericVector adj_se = clone(se);
  for (int i = 0; i < se.size(); ++i) {
    adj_se[i] = std::sqrt(se[i] * se[i] + s_tau2);
  }

  // Boundaries for root-finding
  double min_es = *std::min_element(es.begin(), es.end()) - *std::max_element(se.begin(), se.end());
  double max_es = *std::max_element(es.begin(), es.end()) + *std::max_element(se.begin(), se.end());

  // Generate a random uniform number
  double u = R::runif(0, 1);

  // Define the Edgington combined p-value function as a lambda function
  fntl::dfd f = [&](double h0) {
    return pfct_edge_cpp(NumericVector::create(h0), es, adj_se)[0] - u;
  };

  // Root-finding doesnt break evaluation in case of non-convergence
  fntl::findroot_args args;
  args.action = fntl::error_action::NONE;

  // Find the root using Brent algorithm
  auto out = fntl::findroot_brent(f, min_es, max_es);

  // Return
  return out.root;
}

// Function to generate N samples of mu from the confidence distribution
// [[Rcpp::export]]

NumericVector sample_mu_cpp(NumericVector s_tau2,
                            NumericVector es,
                            NumericVector se) {

  // Initialize a vector for roots
  NumericVector roots(s_tau2.size());

  // Iterate over all samples of tau2
  for (int ii = 0; ii < s_tau2.size(); ++ii) {
    try {
      roots[ii] = sample_one_mu(s_tau2[ii], es, se);
    } catch (...) {
      roots[ii] = NA_REAL;
    }
  }

  // Return
  return roots;
}



// Confidence density of mu (evaluated for one mu)
// [[Rcpp::export]]

double cd_single(NumericVector h0,
                 NumericVector es,
                 NumericVector se,
                 double h) {

  // Function for which derivative is computed
  fntl::dfv f = [&](NumericVector x) {
    return pfct_edge_cpp(x, es, se)[0];
  };

  // Compute and return the numerical derivative
  return fntl::fd_deriv(f, h0, 0, h);
}

// Confidence density of mu (evaluated for multiple mu's)
// [[Rcpp::export]]
NumericVector CD_cpp(NumericVector h0,
                     NumericVector es,
                     NumericVector se,
                     double h = 1e-4) { // same as in numDeriv::grad

  NumericVector dv(h0.size());

  for (int ii = 0; ii < h0.size(); ++ii) {
    dv[ii] = cd_single(NumericVector::create(h0[ii]), es, se, h);
  }

  return dv;

}











