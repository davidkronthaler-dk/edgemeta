// Functions for density of tau2 based on the generalized heterogeneity statistic             

// [[Rcpp::depends(fntl)]]
#include <Rcpp.h>
#include "fntl.h"

#include "metaprediction.h"

using namespace Rcpp;
using namespace std;


// Generalized Q statistic
// [[Rcpp::export]]
double Q_cpp(Rcpp::NumericVector es, 
             Rcpp::NumericVector se, 
             double tau2) {
  
  // Adjusted standard errors
  Rcpp::NumericVector sea(se.size());
  for (int ii = 0; ii < se.size(); ++ii) {
    sea[ii] = std::sqrt(se[ii] * se[ii] + tau2);
  }
  
  double opt = opti_edge(es, sea);
  
  double Q = 0.0;
  
  for (int ll = 0; ll < es.size(); ++ll) {
    Q += std::pow(tau2 + std::pow(se[ll], 2), -1) * std::pow(es[ll] - opt, 2);
  }
  
  return Q;
  
}

// Derivative of generalized Q-statistic
// [[Rcpp::export]]
double dQ_cpp(Rcpp::NumericVector es,
              Rcpp::NumericVector se,
              double tau2,
              double h = 1e-4) {
  
  fntl::dfv f = [=](Rcpp::NumericVector x) {
    return Q_cpp(es, se, x[0]);  
  };
  
  Rcpp::NumericVector x0 = Rcpp::NumericVector::create(tau2); 
  
  return fntl::fd_deriv(f, x0, 0, h);
}

// Density of tau2 by change of variables
// [[Rcpp::export]]
Rcpp::NumericVector ftau2_cpp(Rcpp::NumericVector es, 
                              Rcpp::NumericVector se, 
                              Rcpp::NumericVector tau2) {
  
  Rcpp::NumericVector fd(tau2.size());
  
  for (int ll = 0; ll < tau2.size(); ++ll) {
    double tau = tau2[ll];
    
    if (tau < 0.0) {
      fd[ll] = 0.0;
      continue;
    }
    
    double Q = Q_cpp(es, se, tau);
    double dQ = dQ_cpp(es, se, tau);
    
    fd[ll] = R::dchisq(Q, es.size() - 1, false) * std::abs(dQ);
  }
  
  return fd;
}





























