// Functions for confidence density of tau2 based on the generalized heterogeneity statistic       

// [[Rcpp::depends(fntl)]]
#include <Rcpp.h>
#include "fntl.h"

#include "metaprediction.h"

using namespace Rcpp;
using namespace std;

// F: generalized Q statistic
// [[Rcpp::export]]
double Q_cpp(Rcpp::NumericVector es, 
             Rcpp::NumericVector se, 
             double tau2) {
  
  // Adjusted standard errors
  Rcpp::NumericVector sea(se.size());
  for (int ii = 0; ii < se.size(); ++ii) {
    sea[ii] = (se[ii] * se[ii] + tau2);
  }
  
  // IVWE for computation of Q
  double sumw = 0.0;
  for (int ii = 0; ii < se.size(); ++ii) {
    sumw += 1.0 / sea[ii];
  }
  
  double opt = 0.0;
  for (int ll = 0; ll < es.size(); ++ll) {
    opt += (1 / sea[ll]) * es[ll];
  }
  opt = opt / sumw;
  
  double Q = 0.0;
  
  for (int ll = 0; ll < es.size(); ++ll) {
    Q += std::pow(tau2 + std::pow(se[ll], 2), -1) * std::pow(es[ll] - opt, 2);
  }
  
  return Q;
  
}

// F: numerical derivative of generalized Q-statistic 
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

// F: analytic derivative of generalized Q-statistic for ivwe
// [[Rcpp::export]]
double dQIVWE(NumericVector es,
              NumericVector se,
              double tau2) {
  
  int k = es.size(); // Number of studies
  
  double sw = 0.0;     // sum of weights
  double swe = 0.0;    // sum of weights * effect sizes
  double da = 0.0;     // derivative numerator 
  double db = 0.0;     // derivative denominator
  
  // Compute weights, their derivatives, and derivative of hmu
  std::vector<double> w(k), dw(k);
  for (int i = 0; i < k; ++i) {
    double vi = se[i] * se[i];
    double denom = vi + tau2;
    
    w[i] = 1.0 / denom;
    dw[i] = -1.0 / (denom * denom);
    
    sw += w[i];
    swe += w[i] * es[i];
    
    da += dw[i] * es[i];
    db += dw[i];
  }
  
  double hmu = swe / sw;
  double dhmu = (sw * da - swe * db) / (sw * sw);
  
  // Compute derivative of Q
  double dQ = 0.0;
  for (int i = 0; i < k; ++i) {
    double diff = es[i] - hmu;
    dQ += dw[i] * diff * diff + 2.0 * w[i] * diff * dhmu;
  }
  
  return dQ;
}

// F: density of tau2 by change of variables
// [[Rcpp::export]]
Rcpp::NumericVector ftau2(Rcpp::NumericVector es, 
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
    double dQ = dQIVWE(es, se, tau);
    // double dQ = dQ_cpp(es, se, tau);
    
    fd[ll] = R::dchisq(Q, es.size() - 1, false) * std::abs(dQ);
  }
  
  return fd;
}
