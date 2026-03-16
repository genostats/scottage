#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include "nullObject.h"

template<typename scalar_t>
scalar_t score_test(Rcpp::XPtr<nullObject<scalar_t>> nm, Rcpp::NumericVector X) {
  int m = nm->T.size();
  int n = nm->start.size();
  int nb_cluster = nm->W.cols();

  if(X.size() != n) Rcpp::stop("length(X) is not ok");
  
  X = Rcpp::ifelse(Rcpp::is_na(X), Rcpp::mean(Rcpp::na_omit(X)), X);
  
  std::vector<scalar_t> S1_gamma(m, 0);
  for(int jk = 0; jk < n; jk++) {
    for(int z : nm->L_times[jk]) S1_gamma[z - 0] += X[jk] * nm->exp_beta_Z[jk];  // MINUS 1
  }
 
  std::vector<scalar_t> Ug(n, 0);
  for(int jk : nm->events) {
    int s = nm->S[jk - 0] - 0;                // MINUS 1
    Ug[jk - 0] = X[jk - 0] - S1_gamma[s]/nm->S0[s];  // MINUS 1
  }

  std::vector<scalar_t> EUg(n, 0);
  for(int jk = 0; jk < n; jk++) {
    for(int l : nm->L_events[jk]) {
      int s = nm->S[l - 0] - 0;    // MINUS 1
      EUg[jk] += nm->exp_beta_Z[jk]/nm->S0[s] * (X[jk] - S1_gamma[s]/nm->S0[s]);
    } 
  }

  std::vector<scalar_t> Wg0(n);
  for(int jk = 0; jk < n; jk++) Wg0[jk] = Ug[jk] - EUg[jk];

  Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> Wg = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>::Zero(nb_cluster);
  for(int jk = 0; jk < n; jk++) {
    Wg[ nm->id[jk] - 0 ] += Wg0[jk];    // MINUS 1
  }

  scalar_t Vgg = 0;
  scalar_t score = 0;

  for(scalar_t x : Wg) {
    score += x;
    Vgg += x*x;
  }

  Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> Vkg = (nm->W) * (Wg);
  scalar_t variance = Vgg - (Vkg.transpose() * (nm->Vkk_i) * Vkg)(0,0);

  return(score*score/variance);
}

// [[Rcpp::export]]
float score_test_float(Rcpp::XPtr<nullObject<float>> nm, Rcpp::NumericVector X) {
  return score_test<float>(nm, X);
}

// [[Rcpp::export]]
double score_test_double(Rcpp::XPtr<nullObject<double>> nm, Rcpp::NumericVector X) {
  return score_test<double>(nm, X);
}
