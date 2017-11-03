#include "simol/statphys/potential/Shearing.hpp"

namespace simol
{

  Shearing::Shearing(Input const & input):
    Potential(input),
    amplitude_(input.amplitude()),
    pulsation_(2 * M_PI / input.length())
  {}
  
  DVec Shearing::totalForce(DVec const& position) const
  {
    DVec force = DVec::Zero(dimension_);
    force(1) = nonEqAmplitude() * sin(pulsation_ * position(0));
    return force;
  }

  double Shearing::operator()(DVec const& /*position*/) const
  {
    return 0;
    //return amplitude_ * (1 - cos(pulsation_ * position(0)) * cos(pulsation_ * position(1)) * cos(pulsation_ * position(2)));
  }

  DVec Shearing::gradient(DVec const& /*position*/) const
  {
    return DVec::Constant(dimension_, 0);
  }

  double Shearing::laplacian(DVec const& /*position*/) const
  {
    return 0;
  }

}
