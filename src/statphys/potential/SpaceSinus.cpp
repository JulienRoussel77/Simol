#include "simol/statphys/potential/SpaceSinus.hpp"

namespace simol
{

  SpaceSinus::SpaceSinus(Input const & input):
    Potential(input),
    amplitude_(input.amplitude()),
    pulsation_(2 * M_PI / input.length())
  {}

  double SpaceSinus::operator()(Vector<double> const& position) const
  {
    return amplitude_ * (1 - cos(pulsation_ * position(0)) * cos(pulsation_ * position(1)) * cos(pulsation_ * position(2)));
  }

  Vector<double> SpaceSinus::gradient(Vector<double> const& position) const
  {
    Vector<double> grad(3);
    grad(0) = sin(pulsation_ * position(0)) * cos(pulsation_ * position(1)) * cos(pulsation_ * position(2));
    grad(1) = cos(pulsation_ * position(0)) * sin(pulsation_ * position(1)) * cos(pulsation_ * position(2));
    grad(2) = cos(pulsation_ * position(0)) * cos(pulsation_ * position(1)) * sin(pulsation_ * position(2));
    grad *= amplitude_ * pulsation_;
    return grad;
  }

  double SpaceSinus::laplacian(Vector<double> const& position) const
  {
    return 3 * pow(pulsation_, 2) * value(position);
  }

}
