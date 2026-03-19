#ifndef _scottage_mean_vector_
#define _scottage_mean_vector_

template<typename vec>
double mean_vector(vec & X) {
  double S = 0;  // on fait le calcul en double meme si c'est un vecteur de floats...
  unsigned int n = 0;
  for(auto a : X) {
    if(!std::isnan(a)) {
      S += a;
      n++;
    }
  }
  if(n == 0) throw std::runtime_error("In mean_vector, all values are NaN");
  return S/n;
}


#endif
