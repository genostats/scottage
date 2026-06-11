#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include "nullObject.h"
#include "score_test_AG_TEST.h"

// [[Rcpp::export]]
float score_test_float_TEST(Rcpp::XPtr<nullObject<float>> nm, Rcpp::NumericVector X) {
  return score_test_AG_TEST<float>(nm, X);
}

// [[Rcpp::export]]
double score_test_double_TEST(Rcpp::XPtr<nullObject<double>> nm, Rcpp::NumericVector X) {
  return score_test_AG_TEST<double>(nm, X);
}
