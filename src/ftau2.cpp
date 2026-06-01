// Approximate confidence density of tau2 by the generalized heterogeneity statistic       

// [[Rcpp::depends(fntl)]]
#include <Rcpp.h>
#include "fntl.h"

#include "edgemeta.h"

using namespace Rcpp;
using namespace std;

// Generalized heterogeneity statistic Q(tau2)
// [[Rcpp::export]]
double Q_cpp(Rcpp::NumericVector es, 
             Rcpp::NumericVector se, 
             double tau2) {
  
  // Adjusted standard errors: adj_se = sqrt(se^2 + tau2)
  Rcpp::NumericVector adj_se(se.size());
  for (int ii = 0; ii < se.size(); ++ii) {
    adj_se[ii] = (se[ii] * se[ii] + tau2);
  }
  
  // Classical inverse-variance weights estimator
  double sumw = 0.0;
  for (int ii = 0; ii < se.size(); ++ii) {
    sumw += 1.0 / adj_se[ii];
  }
  double opt = 0.0;
  for (int ll = 0; ll < es.size(); ++ll) {
    opt += (1 / adj_se[ll]) * es[ll];
  }
  opt = opt / sumw;
  
  // Q(tau2)
  double Q = 0.0;
  for (int ll = 0; ll < es.size(); ++ll) {
    Q += std::pow(adj_se[ll], -1) * std::pow(es[ll] - opt, 2);
  }
  
  return Q;
}

// Analytic derivative of Q(tau2)
// [[Rcpp::export]]
double dQIVWE(NumericVector es,
              NumericVector se,
              double tau2) {
  
  int k = es.size();   // Number of studies
  double sw = 0.0;     // sum of weights
  double swe = 0.0;    // sum of weights * effect sizes
  double da = 0.0;     // derivative numerator 
  double db = 0.0;     // derivative denominator
  
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
  
  double dQ = 0.0;
  for (int i = 0; i < k; ++i) {
    double diff = es[i] - hmu;
    dQ += dw[i] * diff * diff + 2.0 * w[i] * diff * dhmu;
  }
  
  return dQ;
}

// Confidence density c(tau2) by change of variables
// [[Rcpp::export]]
Rcpp::NumericVector ftau2(Rcpp::NumericVector es, 
                          Rcpp::NumericVector se, 
                          Rcpp::NumericVector tau2) {
  
  Rcpp::NumericVector fd(tau2.size());
  
  for (int ll = 0; ll < tau2.size(); ++ll) {
    double tau = tau2[ll];
    
    double Q = Q_cpp(es, se, tau);
    double dQ = dQIVWE(es, se, tau);

    fd[ll] = R::dchisq(Q, es.size() - 1, false) * std::abs(dQ);
  }
  
  return fd;
}

