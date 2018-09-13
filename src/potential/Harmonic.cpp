#include "simol/potential/Harmonic.hpp"

namespace simol
{

  Harmonic::Harmonic(Input const & input):
    Potential(input),
    stiffness_(input.potentialStiffness())
  {}

  ///
  ///sampling a gaussian
  double Harmonic::drawLaw(double localBeta, std::shared_ptr<RNG>& rng, int /*type*/) const
  {
    return center_ + rng->scalarGaussian() / sqrt(localBeta*stiffness_);
  }

  double Harmonic::operator()(double distance) const
  { return stiffness_ / 2 * pow(distance - center_, 2); }
  
  double Harmonic::symmetricValue(double position) const
  {
    return value(position);
  }
  
  double Harmonic::skewsymmetricValue(double /*position*/) const
  {
    return 0;
  }

  double Harmonic::scalarGradient(double distance) const
  {
    return stiffness_ * (distance - center_);
    //return DVec::Constant(1, stiffness_ * (distance - sigmaPot_));
  }

  double Harmonic::laplacian(double /*distance*/) const
  {
    return stiffness_;
  }

  /*double Harmonic::drawLaw(double localBeta, std::shared_ptr<RNG>& rng) const
  {
    return center_ + rng->scalarGaussian() / sqrt(localBeta);
  }*/
  
  DVec Harmonic::polyCoeffs() const
  {
    DVec coeffs(DVec::Zero(3));
    
    coeffs[2] = 1./2;
    
    return stiffness_ * coeffs;
  }

}
