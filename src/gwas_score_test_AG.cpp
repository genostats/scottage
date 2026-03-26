#include <Rcpp.h>
#include "snp_filler_with_repeats_bed.h"
#include "gwas_score_test_AG.h"

// using namespace Rcpp;

template<typename filler>
inline List GWAS_score_test_AG_bed(Rcpp::XPtr<nullObject<typename filler::type>> nm, filler S) {
  gwas_score_test_AG<typename filler::type> x(nm, S);
  x.run_tests();
  return S.L;
}


//[[Rcpp::export]]
List GWAS_score_test_AG_bed_float(XPtr<matrix4> pA, NumericVector p, IntegerVector times, int beg, int end, Rcpp::XPtr<nullObject<float>> nm) {
  snp_filler_with_repeats_bed<float> S(pA, p, times, beg, end);
  return GWAS_score_test_AG_bed(nm, S);
}

//[[Rcpp::export]]
List GWAS_score_test_AG_bed_double(XPtr<matrix4> pA, NumericVector p, IntegerVector times, int beg, int end, Rcpp::XPtr<nullObject<double>> nm) {
  snp_filler_with_repeats_bed<double> S(pA, p, times, beg, end);
  return GWAS_score_test_AG_bed(nm, S);
}

