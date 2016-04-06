#include "Sinusoidal.hpp"

namespace simol
{

  Sinusoidal::Sinusoidal(Input const & input):
    Potential(input), 
    amplitude_(input.amplitude()), 
    pulsation_(2*M_PI/input.length())
  {}
  
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
    return amplitude_* (1-cos(pulsation_*position)); 
  }
  
  Vector<double> Sinusoidal::gradient(double position) const
  { 
    return Vector<double>(1, amplitude_*pulsation_*sin(pulsation_*position));
  }
  
  double Sinusoidal::laplacian(double position) const
  {
    return amplitude_ * pow(pulsation_, 2) * cos(pulsation_ * position);
  }
  
}