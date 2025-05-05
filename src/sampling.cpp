// [[Rcpp::depends(fntl)]]
#include <Rcpp.h>
#include "fntl.h"

#include "metaprediction.h"

using namespace Rcpp;
using namespace std;


// Function to generate on sample of mu from the confidence distribution
// [[Rcpp::export]]

double sample_one_mu(double s_tau2,
                     Rcpp::NumericVector es,
                     Rcpp::NumericVector se) {
  
  // Adjust the standard errors based on tau2
  Rcpp::NumericVector adj_se = clone(se);
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
    return pfct_edge_cpp(Rcpp::NumericVector::create(h0), es, adj_se)[0] - u;
  };
  
  // Root-finding doesnt break evaluation in case of non-convergence
  fntl::findroot_args args;
  args.action = fntl::error_action::NONE;
  
  // Find the root using Brent algorithm
  auto out = fntl::findroot_brent(f, min_es, max_es, args);
  
  // Return
  return out.root;
}

// Function to generate B samples of mu from the confidence distribution
// [[Rcpp::export]]

Rcpp::NumericVector sample_mu_cpp(Rcpp::NumericVector s_tau2,
                                  Rcpp::NumericVector es,
                                  Rcpp::NumericVector se) {
  
  // Initialize a vector for roots
  Rcpp::NumericVector roots(s_tau2.size());
  
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


// Function to generate one sample of tau2 from Q(tau2)
// [[Rcpp::export]]
double sample_one_tau2(Rcpp::NumericVector es,
                       Rcpp::NumericVector se,
                       double upper) {
  
  
  // Generate samples of Q(tau2)
  double c = R::rchisq(es.size() - 1);
  
  // Define Q(tau2) - C to solve for the root = tau2
  fntl::dfd f = [&](double tau2) {
    return Q_cpp(es, se, tau2) - c;
  };
  
  // Root-finding does break evaluation in case of non-convergence
  fntl::findroot_args args;
  args.action = fntl::error_action::NONE; 
  
  // Find the root using Brent algorithm
  auto out = fntl::findroot_brent(f, 0, upper, args);
  
  // Return
  return out.root;
}


// Function to generate B samples of tau2 from Q(tau2)
// [[Rcpp::export]]

Rcpp::NumericVector sample_tau2_cpp(int ns,
                                    Rcpp::NumericVector es,
                                    Rcpp::NumericVector se,
                                    double upper) {
  
  // Initialize a vector for roots
  Rcpp::NumericVector roots(ns);
  
  // Iterate for N samples
  for (int ii = 0; ii < ns; ++ii) {
    try {
      roots[ii] = sample_one_tau2(es, se, upper);
    } catch (...) {
      roots[ii] = NA_REAL;  // return NA if root-finding fails
    }
  }
  
  // Return
  return roots;
  
} 


