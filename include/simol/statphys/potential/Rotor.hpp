#ifndef SIMOL_ROTOR_HPP
#define SIMOL_ROTOR_HPP

#include "Potential.hpp"

namespace simol
{

  class Rotor : public Potential
  {
    public:
      Rotor(Input const& input);
      double operator()(double position) const;
      double scalarGradient(double position) const;
      double laplacian(double position) const;
      double drawLaw(double localBeta, std::shared_ptr<RNG>& rng_) const;
  };

}

#endif