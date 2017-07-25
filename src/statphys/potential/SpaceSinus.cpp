#include "simol/statphys/potential/SpaceSinus.hpp"

namespace simol
{

  SpaceSinus::SpaceSinus(Input const & input):
    Potential(input),
    amplitude_(input.amplitude()),
    pulsation_(2 * M_PI / input.length())
  {}
  
  DVec SpaceSinus::totalForce(DVec const& position) const
  {
    //cout << "Potential::totalForce" << endl;
    //cout << nonEqForce_ << " - " << gradient(position) << " = " << nonEqForce_ - gradient(position) << endl;
    DVec force = DVec::Zero(dimension_);
    force(1) = nonEqAmplitude() * sin(pulsation_ * position(0));
    return force;
  }

  double SpaceSinus::operator()(DVec const& /*position*/) const
  {
    return 0;
    //return amplitude_ * (1 - cos(pulsation_ * position(0)) * cos(pulsation_ * position(1)) * cos(pulsation_ * position(2)));
  }

  DVec SpaceSinus::gradient(DVec const& /*position*/) const
  {
    return DVec::Constant(dimension_, 0);
    /*DVec grad(3);
    grad(0) = sin(pulsation_ * position(0)) * cos(pulsation_ * position(1)) * cos(pulsation_ * position(2));
    grad(1) = cos(pulsation_ * position(0)) * sin(pulsation_ * position(1)) * cos(pulsation_ * position(2));
    grad(2) = cos(pulsation_ * position(0)) * cos(pulsation_ * position(1)) * sin(pulsation_ * position(2));
    grad *= amplitude_ * pulsation_;
    return grad;*/
  }

  double SpaceSinus::laplacian(DVec const& /*position*/) const
  {
    return 0;
    //return 3 * pow(pulsation_, 2) * value(position);
  }

}
