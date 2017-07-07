#include "simol/statphys/potential/Harmonic.hpp"

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
    /*coeffs[0] = pow(center_, 2)/2;
    coeffs[1] = - center_;
    coeffs[2] = 1./2;*/
    
    coeffs[2] = 1./2;
    
    return stiffness_ * coeffs;
  }
  
  HarmonicFE::HarmonicFE(Input const & input):
    Harmonic(input),
    dimension_(input.potentialDimension()),
    epsilon_(input.galerkinEpsilon())
  {}
  
  double HarmonicFE::entropicPotential(double distance) const
  {
    return -(dimension_-1) * log( sqrt(pow(distance, 2) + pow(epsilon_, 2))/2 + distance/2);
  }
  
  double HarmonicFE::operator()(double distance) const
  { 
      return Harmonic::operator()(distance) + entropicPotential(distance);
  }
  
  double HarmonicFE::symmetricValue(double position) const
  {
    return Harmonic::operator()(position) + (entropicPotential(position)+entropicPotential(2*center_-position))/2;
  }
  
  double HarmonicFE::skewsymmetricValue(double position) const
  {
    return (entropicPotential(position)-entropicPotential(2*center_-position))/2;
  }

  double HarmonicFE::scalarGradient(double distance) const
  {
    return Harmonic::scalarGradient(distance) - (dimension_-1) / sqrt(pow(distance, 2) + pow(epsilon_, 2));
    //return DVec::Constant(1, stiffness_ * (distance - sigmaPot_));
  }

  double HarmonicFE::laplacian(double distance) const
  {
    return Harmonic::laplacian(distance) + (dimension_-1) * distance * pow(pow(distance, 2) + pow(epsilon_, 2), -1.5);
  }
  
  double HarmonicFE::inverseCoeff() const
  {
    return dimension_-1;
  }

}
