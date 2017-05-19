#include "simol/statphys/potential/Sinusoidal.hpp"

namespace simol
{

  Sinusoidal::Sinusoidal(Input const & input):
    Potential(input),
    amplitude_(input.amplitude()),
    pulsation_(2 * M_PI / input.length())
  {
    domainSize() = 2 * M_PI /pulsation_;
  }
  
  const double& Sinusoidal::parameter1() const
  {
    return amplitude_;
  }

  double Sinusoidal::drawLaw(double localBeta, std::shared_ptr<RNG>& rng) const
  {
    bool reject = true;
    double xdraw, udraw;
    int count = 0;
    while (reject)
    {
      xdraw = -M_PI + 2 * rng->scalarUniform() * M_PI;
      udraw = rng->scalarUniform();

      reject = (udraw > exp(- localBeta * (value(xdraw) + 2)));
      count++;
    }
    return xdraw;
  }

  double Sinusoidal::operator()(double position) const
  {
    return amplitude_ * (1 - cos(pulsation_ * position));
  }

  double Sinusoidal::scalarGradient(double position) const
  {
    return amplitude_ * pulsation_ * sin(pulsation_ * position);
  }

  double Sinusoidal::laplacian(double position) const
  {
    return amplitude_ * pow(pulsation_, 2) * cos(pulsation_ * position);
  }

}
