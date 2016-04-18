#ifndef SIMOL_SINUSOIDAL_HPP
#define SIMOL_SINUSOIDAL_HPP

#include "Potential.hpp"

namespace simol
{

  class Sinusoidal : public Potential
  {
    public:
      Sinusoidal(Input const& input);
      double drawLaw(double localBeta, std::shared_ptr<RNG>& rng) const;
      double operator()(double position) const;
      Vector<double> gradient(double position) const;
      double laplacian(double position) const;

    private:
      double amplitude_;
      double pulsation_;
  };

}

#endif