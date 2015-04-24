#ifndef SIMOL_GAUSSIAN_HPP
#define SIMOL_GAUSSIAN_HPP

#include <random>


namespace simol
{
  template<class ScalarType>
  Vector<ScalarType> normal_distribution(size_t seed,
                                         ScalarType const & mean, 
                                         ScalarType const & standardDeviation)
  {
    std::default_random_engine generator (seed); 
    std::normal_distribution<ScalarType> distribution (mean,standardDeviation);
    return distribution(generator);
  }

}



 #endif
