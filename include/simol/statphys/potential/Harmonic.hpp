#ifndef SIMOL_HARMONIC_HPP
#define SIMOL_HARMONIC_HPP

#include "Potential.hpp"

namespace simol
{

  class Harmonic : public Potential
  {
    public:
      Harmonic(Input const& input);
      virtual double operator()(double position) const;
      virtual double symmetricValue(double position) const;
      virtual double skewsymmetricValue(double position) const;
      virtual DVec gradient(double position) const;
      virtual double laplacian(double position) const;
      virtual double drawLaw(double localBeta, std::shared_ptr<RNG>& rng_, int type) const;
      virtual DVec polyCoeffs() const;
    private:
      double stiffness_;
  };
  
  ///
  ///Represents an effective potential of the form V(q) = stiffness_ / 2 * (distance - center_)^2 + dimension_ * log(q)
  ///where log(q) comes from a "Free Energy" (FE) effect due to the dimension
  class HarmonicFE : public Harmonic
  {
    public:
      HarmonicFE(Input const& input);
      virtual double operator()(double position) const;
      virtual double symmetricValue(double position) const;
      virtual double skewsymmetricValue(double position) const;
      virtual DVec gradient(double position) const;
      virtual double laplacian(double position) const;
      //virtual DVec polynomialCoeffs() const;
      virtual double inverseCoeff() const;
    private:
      double dimension_;
  };

}

#endif