#ifndef SIMOL_FPU_HPP
#define SIMOL_FPU_HPP

#include "Potential.hpp"

namespace simol
{
  ///
  /// Potential V(x) = stiffness x^2 / 2 + alpha x^3 / 3 + beta x^4 / 4
  /// Is convex if and only if alpha^2 < 3 beta
  class FPU : public Potential
  {
    public:
      FPU(Input const& input);
      double operator()(double position) const;
      DVec gradient(double position) const;
      double laplacian(double position) const;
      virtual double shiftToHarmonic() const;
      
      virtual double const& parameter1() const;
      virtual double const& parameter2() const;
      //double drawLaw(double localBeta, std::shared_ptr<RNG>& rng_) const;
    private:
      double stiffness_, alpha_, beta_;
  };

}

#endif