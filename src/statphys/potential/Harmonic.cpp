#include "simol/statphys/potential/Harmonic.hpp"

namespace simol
{
  
  Harmonic::Harmonic(Input const & input):
    Potential(input), 
    stiffness_(input.potentialStiffness()),
    sigmaPot_(input.potentialSigma())
  {}
  
  
  double Harmonic::operator()(double distance) const
  { return stiffness_ / 2* pow(distance - sigmaPot_, 2); }

  Vector<double> Harmonic::gradient(double distance) const
  { 
    return Vector<double>(1, stiffness_ * (distance - sigmaPot_));
  }
  
  double Harmonic::laplacian(double /*distance*/) const
  { 
    return stiffness_;
  }
  
  double Harmonic::drawLaw(double localBeta, std::shared_ptr<RNG>& rng) const
  {
    return rng->scalarGaussian() / sqrt(localBeta);
  }  
  
}
