// [[Rcpp::depends(fntl)]]
#include <Rcpp.h>
#include "fntl.h"

#include "metaprediction.h"

using namespace Rcpp;
using namespace std;

// Optimization of Edgington combined p-value functio
//[[Rcpp::export]]
Rcpp::List opti_num(NumericVector es, 
                    NumericVector se, 
                    bool point = true, 
                    bool ci = true, 
                    double levelci = 0.95) {
  
  double point_est = NA_REAL;
  if (point) {
    point_est = opti_edge(es, se);
  }
  
  double cilower = NA_REAL;
  double ciupper = NA_REAL;
  
  if (ci) {
    double min_es = *std::min_element(es.begin(), es.end()) - *std::max_element(se.begin(), se.end());
    double max_es = *std::max_element(es.begin(), es.end()) + *std::max_element(se.begin(), se.end());
    double alpha = 1.0 - levelci;
    
    fntl::dfd fl = [&](double x) {
      return pfct_edge_cpp(NumericVector::create(x), es, se)[0] - alpha / 2.0;
    };
    
    fntl::dfd fu = [&](double x) {
      return pfct_edge_cpp(NumericVector::create(x), es, se)[0] - (1.0 - alpha / 2.0);
    };
    
    fntl::findroot_args args;
    args.action = fntl::error_action::NONE;
    
    // Find roots
    auto result_lower = fntl::findroot_brent(fl, min_es, max_es, args);
    auto result_upper = fntl::findroot_brent(fu, min_es, max_es, args);
    
    cilower = result_lower.root;
    ciupper = result_upper.root;
  }
  
  return List::create(
    _["point"] = point_est,
    _["cilower"] = cilower,
    _["ciupper"] = ciupper
  );
}
  
  
