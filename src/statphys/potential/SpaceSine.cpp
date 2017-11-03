#include "simol/statphys/potential/SpaceSine.hpp"

namespace simol
{

  SpaceSine::SpaceSine(Input const & input):
    Potential(input),
    amplitude_(input.amplitude()),
    pulsation_(2 * M_PI / input.length()),
    coupling_(input.coupling())
  {
    if (dimension_ != 2)
      throw runtime_error("The SpaceSine potential is only defined in dimension 2!");
    domainSize() = 2 * M_PI /pulsation_;
  }

  double SpaceSine::operator()(DVec const& position) const
  {
    return amplitude_ * (2 - cos(pulsation_ * position(0)) - cos(pulsation_ * position(1))) + coupling_ * (1-cos(pulsation_ * position(0))) * (1-cos(pulsation_*position(1)));
  }

  DVec SpaceSine::gradient(DVec const& position) const
  {
    DVec grad(2);
    grad(0) = amplitude_ * pulsation_ * sin(pulsation_ * position(0)) + coupling_ * pulsation_ * sin(pulsation_ * position(0)) * (1-cos(pulsation_*position(1)));
    grad(1) = amplitude_ * pulsation_ * sin(pulsation_ * position(1)) + coupling_ * pulsation_ * (1-cos(pulsation_ * position(0))) * sin(pulsation_*position(1));

    return grad;
  }

  double SpaceSine::laplacian(DVec const& position) const
  {
    return (amplitude_+coupling_) * pow(pulsation_, 2) * (cos(pulsation_ * position(0)) + cos(pulsation_ * position(1))) - 2 * coupling_  * pow(pulsation_, 2) * cos(pulsation_ * position(0)) * cos(pulsation_*position(1));
  }

}
