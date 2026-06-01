// [[Rcpp::depends(fntl)]]
#include <Rcpp.h>
#include "fntl.h"

#include "edgemeta.h"

using namespace Rcpp;
using namespace std;

// Generate B draws mu* from C(mu*|tau2) ~ U(0,1)
// [[Rcpp::export]]
Rcpp::NumericVector samplemu(Rcpp::NumericVector s_tau2,
                             Rcpp::NumericVector es,
                             Rcpp::NumericVector se,
                             Rcpp::NumericVector w) {
  
  Rcpp::NumericVector roots(s_tau2.size());
  
  for (int ii = 0; ii < s_tau2.size(); ++ii) {
    try {
      roots[ii] = samponemu(s_tau2[ii], es, se, w);
    } catch (...) {
      roots[ii] = NA_REAL;
    }
  }
  
  return roots;
}

// Generate one draw mu* from C(mu*|tau2) ~ U(0,1)
// [[Rcpp::export]]
double samponemu(double s_tau2,
                 Rcpp::NumericVector es,
                 Rcpp::NumericVector se,
                 Rcpp::NumericVector w) {
  
  // Adjust standard errors: adj_se = sqrt(se^2 + tau2)
  Rcpp::NumericVector adj_se = clone(se);
  for (int i = 0; i < se.size(); ++i) {
    adj_se[i] = std::sqrt(se[i] * se[i] + s_tau2);
  }
  
  // Bracketing
  double max_se = *std::max_element(se.begin(), se.end());
  double min_es = *std::min_element(es.begin(), es.end()) - 4 * max_se; 
  double max_es = *std::max_element(es.begin(), es.end()) + 4 * max_se; 
  
  // u ~ U(0,1)
  double u = R::runif(0, 1);
  
  // find root of C(mu*|tau2) = u
  fntl::dfd f = [&](double h0) {
    return pfctedgew(Rcpp::NumericVector::create(h0), es, adj_se, w)[0] - u;
  };
  
  fntl::findroot_args args;
  args.action = fntl::error_action::NONE;
  auto out = fntl::findroot_brent(f, min_es, max_es, args);
  
  return out.root;
}

// Generate B draws tau2* from Q(tau2*) ~ X2_k-1
// [[Rcpp::export]]
Rcpp::NumericVector samptau2(int B,
                             Rcpp::NumericVector es,
                             Rcpp::NumericVector se,
                             double upper) {
  
  Rcpp::NumericVector roots(B);
  
  for (int ii = 0; ii < B; ++ii) {
    try {
      roots[ii] = samponetau2(es, se, upper);
    } catch (...) {
      roots[ii] = NA_REAL;  
    }
  }
  
  return roots;
} 

// Generate one draw tau2* from Q(tau2*) ~ X2_k-1
// [[Rcpp::export]]
double samponetau2(Rcpp::NumericVector es,
                   Rcpp::NumericVector se,
                   double upper) {
  
  // c ~ X2_k-1
  double c = R::rchisq(es.size() - 1);
  
  // find root of Q(tau2*) = c
  fntl::dfd f = [&](double tau2) {
    return Q_cpp(es, se, tau2) - c;
  };
  
  fntl::findroot_args args;
  args.action = fntl::error_action::NONE; 
  auto out = fntl::findroot_brent(f, 0, upper, args);
  
  return out.root;
}


// Generate B draws mu* from C(mu*|tau2_fixed)
//[[Rcpp::export]]
Rcpp::NumericVector samplemusimple(int B, 
                                   double tau2,
                                   Rcpp::NumericVector es,
                                   Rcpp::NumericVector se,
                                   Rcpp::NumericVector w) {
  
  // Adjust standard errors: adj_se = sqrt(se^2 + tau2)
  Rcpp::NumericVector adj_se = clone(se);
  for (int i = 0; i < se.size(); ++i) {
    adj_se[i] = std::sqrt(se[i] * se[i] + tau2);
  }
  
  // Bracketing
  double max_se = *std::max_element(se.begin(), se.end());
  double min_es = *std::min_element(es.begin(), es.end()) - 4 * max_se;
  double max_es = *std::max_element(es.begin(), es.end()) + 4 * max_se;
  
  // Generate n_sample random uniform numbers
  Rcpp::NumericVector u(B);
  for (int ll = 0; ll < B; ++ll){
    u[ll] = R::runif(0, 1);
  }
  
  // find roots of C(mu*|tau2_fixed) = u
  fntl::findroot_args args;
  args.action = fntl::error_action::NONE;
  
  Rcpp::NumericVector smu(B);
  for (int ll = 0; ll < B; ++ll) {
    double u_ll = u[ll];
    
    fntl::dfd f = [&](double h0) {
      return pfctedgew(Rcpp::NumericVector::create(h0), es, adj_se, w)[0] - u_ll;
    };
    
    auto out = fntl::findroot_brent(f, min_es, max_es, args);
    smu[ll] = out.root;
  }
  
  return smu;
}
