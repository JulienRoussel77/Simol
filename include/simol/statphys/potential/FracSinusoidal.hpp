#ifndef SIMOL_FRACSINUSOIDAL_HPP
#define SIMOL_FRACSINUSOIDAL_HPP

#include "Potential.hpp"

namespace simol
{

  class FracSinusoidal : public Potential
  {
    public:
      FracSinusoidal(Input const& input);
      double operator()(double position) const;
      DVec gradient(double position) const;
      virtual double laplacian(double position) const;


    private:
      double amplitude_;
      double pulsation_;
  };

}

#endif