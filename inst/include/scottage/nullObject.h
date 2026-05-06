#include <RcppEigen.h>

#ifndef __scottage_null_object__
#define __scottage_null_object__

template<typename scalar_t>
using MATRIX = Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>;

template<typename scalar_t>
using VECTOR = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;

template<typename scalar_t>
struct nullObject {
  std::vector<scalar_t> start;
  std::vector<scalar_t> stop;
 
  std::vector<unsigned int> id;

  std::vector<scalar_t> T;      // ev.times
  std::vector<unsigned int> S;  // ev.index

  std::vector<scalar_t> S0;
  std::vector<scalar_t> exp_beta_Z;

  MATRIX<scalar_t> W;
  MATRIX<scalar_t> Vkk_i;

  std::vector<unsigned int> events;
  std::vector<std::vector<unsigned int>> L_times;
  std::vector<std::vector<unsigned int>> L_events;
 
  // possibility to add a correlation matrix
  MATRIX<scalar_t> C;

  size_t size() { 
    return start.size(); 
  }
};


#endif

