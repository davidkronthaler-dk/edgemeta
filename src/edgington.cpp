// [[Rcpp::depends(fntl)]]
#include <Rcpp.h>
#include "fntl.h"

#include "metaprediction.h"

using namespace Rcpp;
using namespace std;

// Edgington combined p-value function
// [[Rcpp::export]]
Rcpp::NumericVector pfctedge(Rcpp::NumericVector h0,
                             Rcpp::NumericVector es,
                             Rcpp::NumericVector se) {
  int k = es.size();
  int Lh0 = h0.size();
  
  // Individual study p-value functions
  NumericMatrix ps(Lh0, k);
  for (int j = 0; j < k; ++j) {
    for (int i = 0; i < Lh0; ++i) {
      ps(i, j) = 1.0 - R::pnorm(es[j], h0[i], se[j], 1, 0);
    }
  }
  
  // Edgington combined p-values
  Rcpp::NumericVector pcombined(Lh0);
  for (int i = 0; i < Lh0; ++i) {
    Rcpp::NumericVector row(k);
    
    for (int j = 0; j < k; ++j) {
      row[j] = ps(i, j);
    }
    
    double s = sum(row);
    double pE = 0.0;
    
    if (k < 12) {
      int sumUpper = floor(s);
      for (int ll = 0; ll <= sumUpper; ++ll) {
        pE += std::pow(-1, ll) * R::choose(k, ll) * std::pow(s - ll, k);
      }
      pE /= R::gammafn(k + 1);
    } else {
      double z = std::sqrt(12.0 * k) * (s / k - 0.5);
      pE = R::pnorm(z, 0.0, 1.0, 1, 0);
    }
    
    pcombined[i] = pE;
  }
  
  return pcombined;
}

// One-sided Wald p-value function ("greater" alternative)
// [[Rcpp::export]]
Rcpp::NumericVector p_wald(double x,
                           Rcpp::NumericVector es,
                           Rcpp::NumericVector se) {
  
  Rcpp::NumericVector p = es.size();
  for (int ii = 0; ii < es.size(); ++ii) {
    p[ii] = R::pnorm(es[ii], x, se[ii], false, false);
  }
  
  return p;
}


// Point estimate from Edgington combined p-value function
// [[Rcpp::export]]
double opti_edge(Rcpp::NumericVector es,
                 Rcpp::NumericVector se) {
  
  // Find mean(p) = 0.5
  fntl::dfd f = [&](double x) {
    Rcpp::NumericVector p = p_wald(x, es, se);
    double s = 0.0;
    for (int ii = 0; ii < es.size(); ++ii) {
      s += p[ii];
    }
    s/=es.size();
    return s - 0.5;
  };
  
  // Bracket
  double max_se = *std::max_element(se.begin(), se.end());
  double min_es = *std::min_element(es.begin(), es.end()) - max_se;
  double max_es = *std::max_element(es.begin(), es.end()) + max_se;
  
  // find root
  fntl::findroot_args args;
  args.action = fntl::error_action::NONE;
  auto out = fntl::findroot_brent(f, min_es, max_es, args);
  
  return out.root;
}

// Point estimate and confidence interval from Edgington's p-value function
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
    double max_se = *std::max_element(se.begin(), se.end());
    double min_es = *std::min_element(es.begin(), es.end()) - max_se;
    double max_es = *std::max_element(es.begin(), es.end()) + max_se;
    
    double alpha = 1.0 - levelci;
    
    fntl::dfd fl = [&](double x) {
      return pfctedge(NumericVector::create(x), es, se)[0] - alpha / 2.0;
    };
    
    fntl::dfd fu = [&](double x) {
      return pfctedge(NumericVector::create(x), es, se)[0] - (1.0 - alpha / 2.0);
    };
    
    fntl::findroot_args args;
    args.action = fntl::error_action::NONE;
    
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










