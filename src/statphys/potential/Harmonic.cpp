#include "simol/statphys/potential/Harmonic.hpp"

namespace simol
{

  Harmonic::Harmonic(Input const & input):
    Potential(input),
    stiffness_(input.potentialStiffness()),
    center_(input.potentialCenter())
  {}


  double Harmonic::operator()(double distance) const
  { return stiffness_ / 2 * pow(distance - center_, 2); }

  DVec Harmonic::gradient(double distance) const
  {
    return DVec::Constant(1,1,stiffness_ * (distance - center_));
    //return DVec::Constant(1, stiffness_ * (distance - sigmaPot_));
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
