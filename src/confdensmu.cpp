// [[Rcpp::depends(fntl)]]
#include <Rcpp.h>
#include "fntl.h"

#include "metaprediction.h"

using namespace Rcpp;
using namespace std;

// Confidence density of mu (evaluated for one mu)
// [[Rcpp::export]]

double cd_single(Rcpp::NumericVector h0,
                 Rcpp::NumericVector es,
                 Rcpp::NumericVector se,
                 double h = 1e-4) {
  
  // Function for which derivative is computed
  fntl::dfv f = [&](Rcpp::NumericVector x) {
    return pfctedge(x, es, se)[0];
  };
  
  // Compute and return the numerical derivative
  return fntl::fd_deriv(f, h0, 0, h);
}

// Confidence density of mu (evaluates multiple mu's)
// [[Rcpp::export]]
Rcpp::NumericVector CD_cpp(Rcpp::NumericVector h0,
                           Rcpp::NumericVector es,
                           Rcpp::NumericVector se,
                           double h = 1e-4) { // same as in numDeriv::grad
  
  Rcpp::NumericVector dv(h0.size());
  
  for (int ii = 0; ii < h0.size(); ++ii) {
    dv[ii] = cd_single(Rcpp::NumericVector::create(h0[ii]), es, se, h);
  }
  
  return dv;
  
}
  
  
