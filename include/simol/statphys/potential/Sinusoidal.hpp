#ifndef SIMOL_SINUSOIDAL_HPP
#define SIMOL_SINUSOIDAL_HPP

#include "Potential.hpp"

namespace simol
{

  class Sinusoidal : public Potential
  {
    public:
      Sinusoidal(Input const& input);
      virtual string classname() const {return "Sinusoidal";}
      virtual const double& parameter1() const;
      virtual double drawLaw(double localBeta, std::shared_ptr<RNG>& rng) const;
      virtual double operator()(double position) const;
      virtual double scalarGradient(double position) const;
      virtual double laplacian(double position) const;

    private:
      double amplitude_;
      double pulsation_;
  };

}

#endif