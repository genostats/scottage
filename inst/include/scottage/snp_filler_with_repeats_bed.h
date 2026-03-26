#include <Rcpp.h>
#include "gaston/matrix4.h"
#include "milorGWAS/snp_filler.h"

#ifndef _SCOTTAGE_SNP_FILLER_WITH_REPEATS_BED_
#define _SCOTTAGE_SNP_FILLER_WITH_REPEATS_BED_

using namespace Rcpp;

template<typename scalar_t>
class snp_filler_with_repeats_bed : public snp_filler<scalar_t> {
  public:
  XPtr<matrix4> pA;  
  int ncol, true_ncol;
  NumericVector p; 
  IntegerVector times;
  int beg, end;
  int i;
  snp_filler_with_repeats_bed(XPtr<matrix4> pA_, NumericVector p_, IntegerVector times_, int beg_, int end_)
    : snp_filler<scalar_t>(), pA(pA_), ncol(pA->ncol), true_ncol(pA->true_ncol), 
      p(p_), times(times_), beg(beg_), end(end_), i(beg) { };

  bool snp_fill(scalar_t * SNP) {
    if(i > end) {
      this->monomorphic = true; 
      return false; 
    }
    if( std::isnan(p(i)) || p(i) == 0 || p(i) == 1 ) {
      this->monomorphic = true;
      i++;
      return true;
    }
    unsigned int k = 0;
    uint8_t * snp = pA-> data[i];
    scalar_t mu = 2*p(i);
    for(int ii = 0; ii < true_ncol-1; ii++) {
      uint8_t x = snp[ii];
      for(int ss = 0; ss < 4; ss++) { 
        scalar_t g = ((x&3) != 3)?(x&3):mu;
        for(int z = 0; z < times[4*ii + ss]; z++) {
          SNP[k++] = g;
        }
        x >>= 2;
      }
    }
    { int ii = true_ncol-1;
      uint8_t x = snp[ii];
      for(int ss = 0; ss < 4 && 4*ii+ss < ncol; ss++) {
        scalar_t g = ((x&3) != 3)?(x&3):mu;
        for(int z = 0; z < times[4*ii + ss]; z++) {
          SNP[k++] = g;
        }
        x >>= 2;
      }
    }
    i++;
    this->monomorphic = false;
    return true;
  }
};

#endif
