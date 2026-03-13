#include "nullObject.h"
#include <Rcpp.h>

// templater et faire deux exports, un pour float, l'autre pour double ...

template<typename scalar_t>
Rcpp::XPtr<nullObject<scalar_t>> mk_nullObject(Rcpp::List null_model) {
  nullObject<scalar_t> * pno = new nullObject<scalar_t>;
 
  pno->start = Rcpp::as<std::vector<scalar_t>>(null_model["start"]);
  pno->stop  = Rcpp::as<std::vector<scalar_t>>(null_model["stop"]);
  pno->id    = Rcpp::as<std::vector<unsigned int>>(null_model["id"]);

  pno->T = Rcpp::as<std::vector<scalar_t>>(null_model["ev.times"]);
  pno->S = Rcpp::as<std::vector<unsigned int>>(null_model["ev.index"]);

  pno->S0 = Rcpp::as<std::vector<scalar_t>>(null_model["S0"]);
  pno->exp_beta_Z = Rcpp::as<std::vector<scalar_t>>(null_model["exp.beta.Z"]);

  pno->W = Rcpp::as<Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>>(null_model["W"]);
  pno->Vkk_i = Rcpp::as<Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>>(null_model["inverse.variance"]);

  pno->events = Rcpp::as<std::vector<unsigned int>>(null_model["events"]);
  pno->L_times = Rcpp::as<std::vector<std::vector<unsigned int>>>(null_model["L.times"]);
  pno->L_events = Rcpp::as<std::vector<std::vector<unsigned int>>>(null_model["L.events"]);

  Rcpp::XPtr<nullObject<scalar_t>> res(pno);
  return res; 
}

// [[Rcpp::export]]
Rcpp::XPtr<nullObject<float>> mk_nullObject_float(Rcpp::List null_model) {
  return mk_nullObject<float>(null_model);
}

// [[Rcpp::export]]
Rcpp::XPtr<nullObject<double>> mk_nullObject_double(Rcpp::List null_model) {
  return mk_nullObject<double>(null_model);
}

