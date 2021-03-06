#include "simol/potential/Rotor.hpp"

namespace simol
{

  Rotor::Rotor(Input const & input):
    Potential(input)
  {}


  double Rotor::operator()(double distance) const
  {
    //return pow(distance - 1, 2);
    return 1 - cos(distance);
  }

  double Rotor::scalarGradient(double distance) const
  {
    return sin(distance);
  }

  double Rotor::laplacian(double distance) const
  {
    return cos(distance);
  }

  double Rotor::drawLaw(double localBeta, std::shared_ptr<RNG>& rng) const
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

}
