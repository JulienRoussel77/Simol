#ifndef SIMOL_HARMONIC_HPP
#define SIMOL_HARMONIC_HPP

#include "Potential.hpp"

namespace simol
{

  class Harmonic : public Potential
  {
    public:
      Harmonic(Input const& input);
      double operator()(double position) const;
      DVec gradient(double position) const;
      double laplacian(double position) const;
      double drawLaw(double localBeta, std::shared_ptr<RNG>& rng_) const;
      virtual DVec polynomialCoeffs() const;
    private:
      double stiffness_;
      double center_;
  };

}

#endif