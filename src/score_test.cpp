#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include "nullObject.h"
#include "score_test_AG.h"

// [[Rcpp::export]]
float score_test_float(Rcpp::XPtr<nullObject<float>> nm, Rcpp::NumericVector X) {
  return score_test_AG<float>(nm, X);
}

// [[Rcpp::export]]
double score_test_double(Rcpp::XPtr<nullObject<double>> nm, Rcpp::NumericVector X) {
  return score_test_AG<double>(nm, X);
}
