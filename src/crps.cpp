#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double crps(NumericVector s, NumericVector t) {
  int n = s.size();
  int m = t.size();
  
  double sum_consec_diff = 0.0;
  for (int i = 1; i < n; i++) {
    sum_consec_diff += std::abs(s[i] - s[i - 1]);
  }
  double term2 = sum_consec_diff / (2.0 * (n - 1));
  
  double total_crps = 0.0;
  for (int j = 0; j < m; j++) {
    double sum_abs_diff = 0.0;
    for (int i = 0; i < n; i++) {
      sum_abs_diff += std::abs(s[i] - t[j]);
    }
    double term1 = sum_abs_diff / n;
    
    double crps_j = term1 - term2;
    total_crps += crps_j;
  }
  
  return total_crps / m;
}

