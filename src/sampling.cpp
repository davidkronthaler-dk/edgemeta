// [[Rcpp::depends(fntl)]]
#include <Rcpp.h>
#include "fntl.h"

#include "metaprediction.h"

using namespace Rcpp;
using namespace std;


// F: generate one sample from C(mu|tau2)
// [[Rcpp::export]]
double samponemu(double s_tau2,
                 Rcpp::NumericVector es,
                 Rcpp::NumericVector se) {
  
  // Adjust standard errors
  Rcpp::NumericVector adj_se = clone(se);
  for (int i = 0; i < se.size(); ++i) {
    adj_se[i] = std::sqrt(se[i] * se[i] + s_tau2);
  }
  
  // Boundaries for root-finding
  double max_se = *std::max_element(se.begin(), se.end());
  double min_es = *std::min_element(es.begin(), es.end()) - 4 * max_se; 
  double max_es = *std::max_element(es.begin(), es.end()) + 4 * max_se; 
  
  // Generate a random uniform number
  double u = R::runif(0, 1);
  
  // Edgington combined p-value function as a lambda function and find root
  fntl::dfd f = [&](double h0) {
    return pfctedge(Rcpp::NumericVector::create(h0), es, adj_se)[0] - u;
  };
  fntl::findroot_args args;
  args.action = fntl::error_action::NONE;
  auto out = fntl::findroot_brent(f, min_es, max_es, args);
  
  return out.root;
}

// F: Generate B samples from C(mu|tau2)
// [[Rcpp::export]]

Rcpp::NumericVector samplemu(Rcpp::NumericVector s_tau2,
                             Rcpp::NumericVector es,
                             Rcpp::NumericVector se) {
  
  Rcpp::NumericVector roots(s_tau2.size());
  
  for (int ii = 0; ii < s_tau2.size(); ++ii) {
    try {
      roots[ii] = samponemu(s_tau2[ii], es, se);
    } catch (...) {
      roots[ii] = NA_REAL;
    }
  }
  
  return roots;
}

// F: generate one sample of tau2 from Q(tau2)
// [[Rcpp::export]]

double samponetau2(Rcpp::NumericVector es,
                   Rcpp::NumericVector se,
                   double upper) {
  
  // Generate samples of Q(tau2)
  double c = R::rchisq(es.size() - 1);
  
  // Define Q(tau2) - C to solve for the root = tau2
  fntl::dfd f = [&](double tau2) {
    return Q_cpp(es, se, tau2) - c;
  };
  
  fntl::findroot_args args;
  args.action = fntl::error_action::NONE; 
  auto out = fntl::findroot_brent(f, 0, upper, args);
  
  return out.root;
}

// F: generate B samples of tau2 from Q(tau2)
// [[Rcpp::export]]

Rcpp::NumericVector samptau2(int ns,
                             Rcpp::NumericVector es,
                             Rcpp::NumericVector se,
                             double upper) {
  
  Rcpp::NumericVector roots(ns);
  
  for (int ii = 0; ii < ns; ++ii) {
    try {
      roots[ii] = samponetau2(es, se, upper);
    } catch (...) {
      roots[ii] = NA_REAL;  
    }
  }
  
  return roots;
} 

// F: generate samples of mu from confidence distribution for fixed tau2
//[[Rcpp::export]]

Rcpp::NumericVector samplemusimple(int n_samples, 
                                   double tau2,
                                   Rcpp::NumericVector es,
                                   Rcpp::NumericVector se) {
  
  // Adjust standard errors
  Rcpp::NumericVector adj_se = clone(se);
  for (int i = 0; i < se.size(); ++i) {
    adj_se[i] = std::sqrt(se[i] * se[i] + tau2);
  }
  
  // Boundaries for root-finding
  double max_se = *std::max_element(se.begin(), se.end());
  double min_es = *std::min_element(es.begin(), es.end()) - 4 * max_se;
  double max_es = *std::max_element(es.begin(), es.end()) + 4 * max_se;
  
  // Generate random uniform numbers
  Rcpp::NumericVector u(n_samples);
  for (int ll = 0; ll < n_samples; ++ll){
    u[ll] = R::runif(0, 1);
  }
  
  fntl::findroot_args args;
  args.action = fntl::error_action::NONE;
  
  Rcpp::NumericVector smu(n_samples);
  for (int ll = 0; ll < n_samples; ++ll) {
    double u_ll = u[ll];
    
    fntl::dfd f = [&](double h0) {
      return pfctedge(Rcpp::NumericVector::create(h0), es, adj_se)[0] - u_ll;
    };
    
    auto out = fntl::findroot_brent(f, min_es, max_es, args);
    smu[ll] = out.root;
  }
  
  // Return
  return smu;
}
