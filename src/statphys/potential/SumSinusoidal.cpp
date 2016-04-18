#include "simol/statphys/potential/SumSinusoidal.hpp"

namespace simol
{

  SumSinusoidal::SumSinusoidal(Input const & input):
    Potential(input),
    amplitude_(input.amplitude()),
    pulsation_(2 * M_PI / input.length())
  {}

  double SumSinusoidal::operator()(double position) const
  {
    return amplitude_ * (std::sin(pulsation_ * position)
                         + cos(2 * pulsation_ * position)
                         + cos(3 * pulsation_ * position) / 3);
  }

  Vector<double> SumSinusoidal::gradient(double position) const
  {
    Vector<double> deriv(1);
    deriv(0) = amplitude_ * pulsation_ * (cos(pulsation_ * position)
                                          - 2 * sin(2 * pulsation_ * position)
                                          - sin(3 * pulsation_ * position));
    return deriv;
  }


  double SumSinusoidal::laplacian(double position) const
  {
    return amplitude_ * pow(pulsation_, 2) * (-sin(pulsation_ * position)
           - 4 * cos(2 * pulsation_ * position)
           - 3 * cos(3 * pulsation_ * position));
  }

}
