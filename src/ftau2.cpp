// Functions for density of tau2 based on the generalized heterogeneity statistic             

// [[Rcpp::depends(fntl)]]
#include <Rcpp.h>
#include "fntl.h"
using namespace Rcpp;


// Generalized Q statistic
// [[Rcpp::export]]
double Q_cpp(NumericVector es, 
             NumericVector se, 
             double tau2, 
             double mu_hat) {
  
  double Q = 0.0;
  
  for (int ll = 0; ll < es.size(); ++ll) {
    Q += std::pow(tau2 + std::pow(se[ll], 2), -1) * std::pow(es[ll] - mu_hat, 2);
  }
  
  return Q;
  
}

// Derivative of generalized Q statistic
// [[Rcpp::export]]
double dQ_cpp(NumericVector es, 
              NumericVector se, 
              double tau2, 
              double mu_hat){
  
  double dQ = 0.0;
  
  for (int ll = 0; ll < es.size(); ++ll) {
    dQ += -1 * std::pow(tau2 + std::pow(se[ll], 2), -2) * std::pow(es[ll] - mu_hat, 2);
  }
  
  return dQ;
  
}

// Density of tau2 by change of variables
// [[Rcpp::export]]

NumericVector ftau2_cpp(NumericVector es, 
                        NumericVector se, 
                        NumericVector tau2, 
                        double mu_hat){
  
  NumericVector fd(tau2.size());
  
  
  for (int ll = 0; ll < tau2.size(); ++ll) {
    
    if (tau2[ll] >= 0) {
      
      fd[ll] = R::dchisq(Q_cpp(es, se, tau2[ll], mu_hat), es.size() - 1, false) *
        std::abs(dQ_cpp(es, se, tau2[ll], mu_hat));
      
    } else {
      
      fd[ll] = 0.0;
      
    }
  }
  
  return fd;
  
}


// Function to generate one sample of tau2 from Q(tau2)
// [[Rcpp::export]]


double sample_one_tau2(NumericVector es,
                       NumericVector se,
                       double mu_hat,
                       double upper) {
  
  
  // Generate samples of Q(tau2)
  double c = R::rchisq(es.size() - 1);
  
  // Define Q(tau2) - C to solve for the root = tau2
  fntl::dfd f = [&](double tau2) {
    return Q_cpp(es, se, tau2, mu_hat) - c;
  };
  
  // Root-finding does break evaluation in case of non-convergence
  fntl::findroot_args args;
  args.action = fntl::error_action::NONE; 
  
  // Find the root using Brent algorithm
  auto out = fntl::findroot_brent(f, 0, upper);
  
  // Return
  return out.root;
}


// Function to generate N samples of tau2 from Q(tau2)
// [[Rcpp::export]]

NumericVector sample_tau2_cpp(int ns,
                              NumericVector es,
                              NumericVector se,
                              double mu_hat,
                              double upper) {
  
  // Initialize a vector for roots
  NumericVector roots(ns);
  
  // Iterate for N samples
  for (int ii = 0; ii < ns; ++ii) {
    try {
      roots[ii] = sample_one_tau2(es, se, mu_hat, upper);
    } catch (...) {
      roots[ii] = NA_REAL;  // return NA if root-finding fails
    }
  }
  
  // Return
  return roots;
  
} 





  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  