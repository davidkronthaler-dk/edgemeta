#include <Rcpp.h>
using namespace Rcpp;

// Efficient Monte Carlo approximation to compute univariate CRPS
// [[Rcpp::export]]
double crpsCPP(NumericVector s, NumericVector t) {
  int n = s.size();
  int m = t.size();
  
  double sum_consec_diff = 0.0;
  for (int i = 1; i < n; i++) {
    sum_consec_diff += std::abs(s[i] - s[i - 1]);
  }
  double trm2 = sum_consec_diff / (2.0 * (n - 1));
  
  double sum_crps = 0.0;
  for (int j = 0; j < m; j++) {
    double sum_abs_diff = 0.0;
    for (int i = 0; i < n; i++) {
      sum_abs_diff += std::abs(s[i] - t[j]);
    }
    double trm1 = sum_abs_diff / n;
    
    double crps_j = trm1 - trm2;
    sum_crps += crps_j;
  }
  
  return sum_crps / m;
}

