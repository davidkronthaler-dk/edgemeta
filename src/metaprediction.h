#ifndef MY_FUNCTIONS_H
#define MY_FUNCTIONS_H

#include <Rcpp.h>
using namespace Rcpp;

// Declared functions
Rcpp::NumericVector pfctedge(Rcpp::NumericVector h0, Rcpp::NumericVector es, 
                                  Rcpp::NumericVector se);

Rcpp::NumericVector p_wald(double x, Rcpp::NumericVector es, Rcpp::NumericVector se);

double opti_edge(Rcpp::NumericVector es, Rcpp::NumericVector se);

double cd_single(Rcpp::NumericVector h0, Rcpp::NumericVector es, 
                 Rcpp::NumericVector se, double h);

Rcpp::NumericVector CD_cpp(Rcpp::NumericVector h0, Rcpp::NumericVector es, 
                           Rcpp::NumericVector se, double h);

double samponemu(double s_tau2, Rcpp::NumericVector es, Rcpp::NumericVector se);

Rcpp::NumericVector samplemu(Rcpp::NumericVector s_tau2, Rcpp::NumericVector es, 
                            Rcpp::NumericVector se);

double samponetau2(Rcpp::NumericVector es, Rcpp::NumericVector se, double upper);

Rcpp::NumericVector samptau2(int ns, Rcpp::NumericVector es, 
                                    Rcpp::NumericVector se, double upper);

double Q_cpp(Rcpp::NumericVector es, Rcpp::NumericVector se, double tau2);

double dQ_cpp(Rcpp::NumericVector es, Rcpp::NumericVector se, double tau2, 
              double h) ;

double dQIVWE(NumericVector es, NumericVector se, double tau2);

Rcpp::NumericVector ftau2(Rcpp::NumericVector es, Rcpp::NumericVector se, 
                              Rcpp::NumericVector tau2);

Rcpp::List opti_num(NumericVector es, NumericVector se, bool point, 
                    bool ci, double levelci);

#endif