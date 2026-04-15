// [[Rcpp::depends(fntl)]]
#include <Rcpp.h>
#include "fntl.h"

#include "metaprediction.h"

using namespace Rcpp;
using namespace std;

// Normalizing constant: int c(tau2) dtau2
// [[Rcpp::export]]
double nor_ctau2(NumericVector es, NumericVector se, double utau2) {
  
  fntl::integrate_args args;
  args.subdivisions = 1000L; 
  
  fntl::dfd f = [&](double tau2_scalar) -> double {
    Rcpp::NumericVector tau2_vec(1);
    tau2_vec[0] = tau2_scalar;
    
    Rcpp::NumericVector val = ftau2(es, se, tau2_vec);
    return val[0];
  };
  
  fntl::integrate_result out = fntl::integrate(f, 0, utau2, args);
  return out.value; 
}

// Joint confidence density: c(mu | tau2) * c(tau2)
// [[Rcpp::export]]
double jointCD(Rcpp::NumericVector mu, Rcpp::NumericVector tau2, 
               Rcpp::NumericVector es, Rcpp::NumericVector se,
               double C) {
  int k = se.size();
  Rcpp::NumericVector adj_se(k);
  for (int i = 0; i < k; ++i) {
    adj_se[i] = std::sqrt(se[i] * se[i] + tau2[0]);
  }
  
  Rcpp::NumericVector cdmu = CD_cpp(mu, es, adj_se, 1e-4);
  Rcpp::NumericVector cdtau2 = ftau2(es, se, tau2) / C; // divide by normalizing constant    
  
  return cdmu[0] * cdtau2[0];
}

// Marginal confidence density: c(mu) = int c(mu |tau2) * c(tau2) dtau2
// [[Rcpp::export]]
double marCDsingle(double mu, Rcpp::NumericVector es, Rcpp::NumericVector se,
                   double utau2) {
  
  double C = nor_ctau2(es, se, utau2);
  double ltau2 = 0;
  
  fntl::integrate_args args;
  args.subdivisions = 1000L;
  
  Rcpp::NumericVector mu_vec(1);
  mu_vec[0] = mu;
  
  fntl::dfd f = [&](double tau2_scalar) -> double {
    Rcpp::NumericVector tau2_vec(1);
    tau2_vec[0] = tau2_scalar;
    return jointCD(mu_vec, tau2_vec, es, se, C);
  };
  
  fntl::integrate_result out = fntl::integrate(f, ltau2, utau2, args);
  
  return out.value;
}

// [[Rcpp::export]]
Rcpp::NumericVector marCD(Rcpp::NumericVector mu, Rcpp::NumericVector es,
                          Rcpp::NumericVector se, double utau2) {
  
  Rcpp::NumericVector inter(mu.size());
  for (int ii = 0; ii < mu.size(); ++ii) {
    inter[ii] = marCDsingle(mu[ii], es, se, utau2);
  }
  
  return inter;
}

// Find confidence intervals from c(mu)
// [[Rcpp::export]]
NumericMatrix reff(NumericVector es,
                   NumericVector se,
                   double utau2,
                   double grid_step = 1e-2) {
  
  // Grid
  double max_se = *std::max_element(se.begin(), se.end());
  double lbb = *std::min_element(es.begin(), es.end()) - 4 * max_se;
  double ubb = *std::max_element(es.begin(), es.end()) + 4 * max_se;

  int ngrid = std::ceil((ubb - lbb) / grid_step) + 1;
  NumericVector xgrid(ngrid);
  
  for (int i = 0; i < ngrid; ++i) {
    xgrid[i] = lbb + i * grid_step;
  }  
  
  // Marginal confidence density mu
  NumericVector ff(ngrid);
  for (int i = 0; i < ngrid; ++i) {
    ff[i] = marCDsingle(xgrid[i], es, se, utau2);
  }
  
  // Normalize marginal density (trapezoidal rule)
  double area = 0.0;
  for (int i = 1; i < ngrid; ++i) {
    area += 0.5 * (ff[i] + ff[i - 1]) * grid_step;
  }
  for (int i = 0; i < ngrid; ++i) {
    ff[i] /= area;
  }

  //CDF by numerical integration
  NumericVector cdfn(ngrid);
  cdfn[0] = ff[0] * grid_step;
  
  for (int i = 1; i < ngrid; ++i) {
    // cdfn[i] = cdfn[i - 1] + ff[i] * grid_step; // left Riemanns sum
    cdfn[i] = cdfn[i - 1] + 0.5 * (ff[i] + ff[i - 1]) * grid_step; // trapezoidal rule
    
  }
  NumericVector cdf = cdfn / cdfn[ngrid - 1];  // Normalize to [0,1]
  
  NumericMatrix mat(ngrid, 2);
  for (int i = 0; i < ngrid; ++i) {
    mat(i, 0) = xgrid[i];
    mat(i, 1) = cdf[i];
  }
  
  return mat;
  
}
