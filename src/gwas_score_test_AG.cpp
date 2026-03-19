#include <Rcpp.h>
#include "milorGWAS/snp_filler_bed.h"
#include "milorGWAS/snp_filler_001_bed.h"
#include "milorGWAS/snp_filler_011_bed.h"
#include "gwas_score_test_AG.h"

// using namespace Rcpp;

template<typename filler>
inline List GWAS_score_test_AG_bed(Rcpp::XPtr<nullObject<typename filler::type>> nm, filler S) {
  gwas_score_test_AG<typename filler::type> x(nm, S);
  x.run_tests();
  return S.L;
}


//[[Rcpp::export]]
List GWAS_score_test_AG_bed_float(XPtr<matrix4> pA, NumericVector p, int beg, int end, Rcpp::XPtr<nullObject<float>> nm, std::string coding) {
  if(coding == "012") {
    snp_filler_additive_bed<float> S(pA, p, beg, end);
    return GWAS_score_test_AG_bed(nm, S);
  } else if(coding == "011") {
    snp_filler_011_bed<float> S(pA, p, beg, end);
    return GWAS_score_test_AG_bed(nm, S);
  } else if(coding == "001") {
    snp_filler_001_bed<float> S(pA, p, beg, end);
    return GWAS_score_test_AG_bed(nm, S);
  } else {
    stop("Unknown coding value");
  }
}

