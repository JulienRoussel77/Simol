#ifndef SIMOL_FPU_HPP
#define SIMOL_FPU_HPP

#include "Potential.hpp"

namespace simol
{
  ///
  /// Potential V(x) = stiffness x^2 / 2 + alpha x^3 / 3 + beta x^4 / 4
  /// Is convex if and only if alpha^2 < 3 beta
  /// parameter1_ : alpha
  /// parameter2_ : beta
  class FPU : public Potential
  {
    public:
      FPU(Input const& input);
      virtual string classname() const {return "FPU";}
      void computeHarmonic();
      double operator()(double position) const;
      double scalarGradient(double position) const;
      double laplacian(double position) const;
      virtual double shiftToHarmonic() const;
      
      virtual double const& parameter1() const;
      virtual double const& parameter2() const;
      virtual double harmonicForce(double dist) const;
      virtual double harmonicStiffness() const;
      virtual double harmonicEquilibrium() const;
      virtual double harmonicFrequency() const;
      //double drawLaw(double localBeta, std::shared_ptr<RNG>& rng_) const;
    private:
      double stiffness_, alpha_, beta_;
      double qRepartitionFct_, m1_, m2_, m3_, m4_;
      double harmonicEquilibrium_, harmonicStiffness_;
  };

}

#endif