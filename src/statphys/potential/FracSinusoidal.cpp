#include "simol/statphys/potential/FracSinusoidal.hpp"

namespace simol
{

  FracSinusoidal::FracSinusoidal(Input const & input):
    Potential(input),
    amplitude_(input.amplitude()),
    pulsation_(2 * M_PI / input.length())
  {}

  double FracSinusoidal::operator()(double position) const
  {
    return amplitude_ * cos(2 * M_PI * position) / (2 + sin(2 * M_PI * position));
  }

  DVec FracSinusoidal::gradient(double position) const
  {
    return DVec::Constant(1,1, -amplitude_ * (4 * M_PI * sin(2 * M_PI * position) + 2 * M_PI) / pow(sin(2 * M_PI * position) + 2, 2));
  }


  double FracSinusoidal::laplacian(double position) const
  {
    return -amplitude_ * 32 * pow(M_PI, 2) * pow(sin(M_PI / 4 - M_PI * position), 3) * sin(M_PI * position + M_PI / 4)
           / pow(sin(2 * M_PI * position) + 2, 3);
  }

}
