#include "nullObject.h"
#include <Rcpp.h>

// templater et faire deux exports, un pour float, l'autre pour double ...

// [[Rcpp::export]]
Rcpp::List get_nullObject(Rcpp::XPtr<nullObject<float>> nm) {
  Rcpp::List L;
  L["start"] = Rcpp::wrap(nm->start);
  L["stop"] = Rcpp::wrap(nm->stop);
  L["id"] = Rcpp::wrap(nm->id);
  L["ev.times"] = Rcpp::wrap(nm->T);
  L["ev.index"] = Rcpp::wrap(nm->S);
  L["S0"] = Rcpp::wrap(nm->S0);
  L["exp.beta.Z"] = Rcpp::wrap(nm->exp_beta_Z);
  L["W"] = Rcpp::wrap(nm->W);
  L["inverse.variance"] = Rcpp::wrap(nm->Vkk_i);
  L["events"] = Rcpp::wrap(nm->events);
  L["L.times"] = Rcpp::wrap(nm->L_times);
  L["L.events"] = Rcpp::wrap(nm->L_events);

  return L; 
}
