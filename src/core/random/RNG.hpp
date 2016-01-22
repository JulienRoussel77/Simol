#ifndef SIMOL_GAUSSIAN_HPP
#define SIMOL_GAUSSIAN_HPP

#include <random>
#include "linalg/Vector.hpp"


namespace simol
{
  class RNG
  {
    int dimension_;
    size_t seed_;
    std::mt19937_64 generator_; 
    std::normal_distribution<double> distribution_;
  public:
    RNG(size_t const seed, int dimension):dimension_(dimension), seed_(seed), generator_(seed), distribution_(0,1){}
    
    //Vector<double> gaussian(double const & mean, double const & standardDeviation)
    dvec gaussian()
    {
      dvec g(dimension_);
      for (int i=0; i<dimension_; i++)
	g(i) = distribution_(generator_);
      return g;
    }
  };

}



 #endif
