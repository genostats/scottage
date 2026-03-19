#include <Rcpp.h>
#include "milorGWAS/snp_filler.h"
#include "nullObject.h"
#include "score_test_AG.h"

#ifndef _scottage_gwas_score_test_AG_
#define _scottage_gwas_score_test_AG_

template<typename scalar_t>
class gwas_score_test_AG {
  nullObject<scalar_t> * nm;
  snp_filler<scalar_t> & S;
  std::vector<scalar_t> SNP;

  public:
  gwas_score_test_AG(nullObject<scalar_t> * nm_, snp_filler<scalar_t> & S_) : nm(nm_), S(S_), SNP(nm->size()) {}

  void run_tests() {
    std::vector<scalar_t> score;
    
    while( S.snp_fill( &SNP[0] ) ) {
      if( S.current_snp_monomorphic() ) {
        score.push_back(NAN);
        continue;
      }
      scalar_t s = score_test_AG(nm, SNP);
      score.push_back(s);
    }
    S.L["score"] = Rcpp::wrap(score);
  }

};

#endif
