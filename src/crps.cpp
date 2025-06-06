#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
using namespace Rcpp;

// Function to compute Continuous Ranked Probability Score (CRPS)
// [[Rcpp::export]]
double crps(const std::vector<double>& s, const std::vector<double>& tn) {
  int n = s.size();
  
  // Sort the predictive samples
  std::vector<double> pd_sorted = s;
  std::sort(pd_sorted.begin(), pd_sorted.end());
  
  // Compute cumulative weights: 0, 1, ..., n-1
  std::vector<double> cum_weights(n);
  for (int i = 0; i < n; ++i) {
    cum_weights[i] = i;
  }
  
  // Efficient computation of term2
  // This represents the average pairwise absolute difference between forecast samples
  double term2_a = 0.0;
  double term2_b = 0.0; 
  
  for (int i = 0; i < n; ++i) {
    term2_a += cum_weights[i] * pd_sorted[i];
    term2_b += cum_weights[n - 1 - i] * pd_sorted[i];
  }
  
  // Final term2 calculation (equivalent to sum(|s_i - s_j|) / n^2)
  double term2 = (2.0 / (n * n)) * (term2_a - term2_b);
  
  // Compute term1 for each true value in tn
  // term1 is the average absolute error between forecast samples and each true value
  double crps_sum = 0.0;
  for (double t : tn) {
    double term1 = 0.0;
    for (double sample : s) {
      term1 += std::abs(sample - t);
    }
    term1 /= n;
    
    // Add CRPS for this observation: term1 - 0.5 * term2
    crps_sum += (term1 - 0.5 * term2);
  }
  
  // return average CRPS across all observations
  return crps_sum / tn.size();
}




