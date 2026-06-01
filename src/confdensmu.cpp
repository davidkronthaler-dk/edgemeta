// [[Rcpp::depends(fntl)]]
#include <Rcpp.h>
#include "fntl.h"

#include "edgemeta.h"

using namespace Rcpp;
using namespace std;

// Confidence density of mu; numeric derivative of Edgington's p-value function
// [[Rcpp::export]]
Rcpp::NumericVector CD_cpp(Rcpp::NumericVector h0,
                           Rcpp::NumericVector es,
                           Rcpp::NumericVector se,
                           Rcpp::NumericVector w,
                           double h = 1e-4) {
  
  fntl::dfv f = [&](Rcpp::NumericVector x) {
    return pfctedgew(x, es, se, w)[0];
  };
  
  Rcpp::NumericVector dv(h0.size());
  
  for (int ii = 0; ii < h0.size(); ++ii) {
    dv[ii] = fntl::fd_deriv(f, Rcpp::NumericVector::create(h0[ii]), 0, h);
  }
  
  return dv;
}
