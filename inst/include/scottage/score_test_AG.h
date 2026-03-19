#include <iostream>
#include "nullObject.h"
#include "mean_vector.h"

#ifndef _scottage_core_test_AG_
#define _scottage_core_test_AG_

template<typename scalar_t, typename vector>
scalar_t score_test_AG(nullObject<scalar_t> * nm, vector & X, scalar_t mean_X = std::numeric_limits<scalar_t>::quiet_NaN()) {
  int m = nm->T.size();
  int n = nm->start.size();
  int nb_cluster = nm->W.cols();

  if(X.size() != n) throw std::runtime_error("length(X) is not ok");
  
  std::vector<scalar_t> S1_gamma(m, 0);
  for(int jk = 0; jk < n; jk++) {
    // dans la première boucle on se préoccupe de calculer mean_X si besoin
    scalar_t Xjk;
    if(std::isnan(X[jk])) {
      if(std::isnan(mean_X)) mean_X = mean_vector(X);
      Xjk = mean_X;
    } else {
      Xjk = X[jk];   // ON POURRAIT ACCEPTER DE MODIFIER X "in place" ... ça gagnerait du temps dans la suite (plus besoin de tester)
    }
    for(int z : nm->L_times[jk]) S1_gamma[z - 0] += Xjk * nm->exp_beta_Z[jk];  // MINUS 1
  }
 
  std::vector<scalar_t> Ug(n, 0);
  for(int jk : nm->events) {
    scalar_t Xjk = std::isnan(X[jk])?mean_X:X[jk];
    int s = nm->S[jk];               
    Ug[jk - 0] = Xjk - S1_gamma[s]/nm->S0[s];  // MINUS 1
  }

  std::vector<scalar_t> EUg(n, 0);
  for(int jk = 0; jk < n; jk++) {
    scalar_t Xjk = std::isnan(X[jk])?mean_X:X[jk];
    for(int l : nm->L_events[jk]) {
      int s = nm->S[l];    
      EUg[jk] += nm->exp_beta_Z[jk]/nm->S0[s] * (Xjk - S1_gamma[s]/nm->S0[s]);
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


#endif
