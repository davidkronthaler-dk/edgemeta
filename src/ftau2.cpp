// Functions for density of tau2 based on the generalized heterogeneity statistic             

// [[Rcpp::depends(fntl)]]
#include <Rcpp.h>
#include "fntl.h"
using namespace Rcpp;

// One-sided Wald p-value function
// [[Rcpp::export]]
NumericVector p_wald(double x,
                     NumericVector es,
                     NumericVector se) {
  
  NumericVector p = es.size();
  
  for (int ii = 0; ii < es.size(); ++ii) {
    p[ii] = R::pnorm(es[ii], x, se[ii], false, false);
  }
  
  return p;
}

// Optimize Edgington p-value function
// [[Rcpp::export]]
double opti_edge(NumericVector es,
                 NumericVector se) {
  
  // Find mean(p) = 0.5
  fntl::dfd f = [&](double x) {
    
    // One-sided p-values
    NumericVector p = p_wald(x, es, se);
    
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
  auto out = fntl::findroot_brent(f, min_es, max_es);
  
  // Return
  return out.root;
  
}


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
                        NumericVector tau2) {
  
  NumericVector fd(tau2.size());
  
  for (int ll = 0; ll < tau2.size(); ++ll) {
    double tau = tau2[ll];
    
    if (tau < 0.0) {
      fd[ll] = 0.0;
      continue;
    }
    
    // Adjusted standard errors
    NumericVector sea(se.size());
    for (int ii = 0; ii < se.size(); ++ii) {
      sea[ii] = std::sqrt(se[ii] * se[ii] + tau); 
    }
    
    double mu_hat = opti_edge(es, sea);
    double Q = Q_cpp(es, se, tau, mu_hat);
    double dQ = dQ_cpp(es, se, tau, mu_hat);
    
    fd[ll] = R::dchisq(Q, es.size() - 1, false) * std::abs(dQ);
  }
  
  return fd;
}

// NumericVector ftau2_cpp(NumericVector es, 
//                         NumericVector se, 
//                         NumericVector tau2,
//                         double mu_hat){
//   
//   
//   NumericVector fd(tau2.size());
//   
//   for (int ll = 0; ll < tau2.size(); ++ll) {
//     
//     if (tau2[ll] >= 0) {
//       
//       fd[ll] = R::dchisq(Q_cpp(es, se, tau2[ll], mu_hat), es.size() - 1, false) *
//         std::abs(dQ_cpp(es, se, tau2[ll], mu_hat));
//       
//     } else {
//       
//       fd[ll] = 0.0;
//       
//     }
//   }
//   
//   return fd;
//   
// }

// Function to generate one sample of tau2 from Q(tau2)
// [[Rcpp::export]]
double sample_one_tau2(NumericVector es,
                       NumericVector se,
                       double upper) {
  
  
  // Generate samples of Q(tau2)
  double c = R::rchisq(es.size() - 1);
  
  // Define Q(tau2) - C to solve for the root = tau2
  fntl::dfd f = [&](double tau2) {
    
    // Adjusted standard errors
    NumericVector sea(se.size());
    for (int ii = 0; ii < se.size(); ++ii) {
      sea[ii] = std::sqrt(se[ii] * se[ii] + tau2); 
    }
    
    double mu_hat = opti_edge(es, sea);
    
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
                              double upper) {
  
  // Initialize a vector for roots
  NumericVector roots(ns);
  
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





  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  