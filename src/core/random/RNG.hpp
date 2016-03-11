#ifndef SIMOL_GAUSSIAN_HPP
#define SIMOL_GAUSSIAN_HPP

#include <random>
#include "core/linalg/Vector.hpp"


namespace simol
{
  class RNG
  {
    int dimension_;
    size_t seed_;
    std::mt19937_64 generator_; 
    std::normal_distribution<double> normalDistribution_;
		std::uniform_real_distribution<double> uniformDistribution_;
		//std::poisson_distribution<int> poissonDistribution_;
		//std::geometric_distribution<int> geometricDistribution_;
		std::exponential_distribution<double> exponentialDistribution_;
  public:
    RNG(size_t const seed, int dimension):
			dimension_(dimension), 
			seed_(seed), 
			generator_(seed), 
			normalDistribution_(0,1),
			uniformDistribution_(0.,1.),
			exponentialDistribution_(1)
		{}
    
    //Vector<double> gaussian(double const & mean, double const & standardDeviation)
    dvec gaussian()
    {
      dvec g(dimension_);
      for (int i=0; i<dimension_; i++)
				g(i) = normalDistribution_(generator_);
      return g;
    }
    
    double scalarGaussian()
    {
      return normalDistribution_(generator_);
    }
    
    dvec uniform()
    {
      dvec u(dimension_);
      for (int i=0; i<dimension_; i++)
				u(i) = uniformDistribution_(generator_);
      return u;
    }
    
    double scalarUniform()
		{
			return uniformDistribution_(generator_);
		}
		
		/*int scalarPoisson()
		{
			return poissonDistribution_(generator_);
		}
		
		int scalarGeometric()
		{
			return geometricDistribution_(generator_);
		}*/
		
    double scalarExponential()
		{
			return exponentialDistribution_(generator_);
		}
		
  };

}



 #endif
