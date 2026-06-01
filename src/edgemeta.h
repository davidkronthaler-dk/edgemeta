#ifndef MY_FUNCTIONS_H
#define MY_FUNCTIONS_H

#include <Rcpp.h>

using namespace Rcpp;

Rcpp::NumericVector pfctedgew(
    Rcpp::NumericVector h0,
    Rcpp::NumericVector es,
    Rcpp::NumericVector se,
    Rcpp::NumericVector w
);

double cd_single(
    Rcpp::NumericVector h0,
    Rcpp::NumericVector es,
    Rcpp::NumericVector se,
    double h
);

Rcpp::NumericVector CD_cpp(
    Rcpp::NumericVector h0,
    Rcpp::NumericVector es,
    Rcpp::NumericVector se,
    Rcpp::NumericVector w,
    double h
);

double samponemu(
    double s_tau2,
    Rcpp::NumericVector es,
    Rcpp::NumericVector se,
    Rcpp::NumericVector w
);

Rcpp::NumericVector samplemu(
    Rcpp::NumericVector s_tau2,
    Rcpp::NumericVector es,
    Rcpp::NumericVector se,
    Rcpp::NumericVector w
);

double samponetau2(
    Rcpp::NumericVector es,
    Rcpp::NumericVector se,
    double upper
);

Rcpp::NumericVector samptau2(
    int ns,
    Rcpp::NumericVector es,
    Rcpp::NumericVector se,
    double upper
);

double Q_cpp(
    Rcpp::NumericVector es,
    Rcpp::NumericVector se,
    double tau2
);

double dQIVWE(
    Rcpp::NumericVector es,
    Rcpp::NumericVector se,
    double tau2
);

Rcpp::NumericVector ftau2(
    Rcpp::NumericVector es,
    Rcpp::NumericVector se,
    Rcpp::NumericVector tau2
);

#endif