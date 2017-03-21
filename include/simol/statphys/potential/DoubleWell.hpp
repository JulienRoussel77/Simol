#ifndef SIMOL_DOUBLEWELL_HPP
#define SIMOL_DOUBLEWELL_HPP

#include "Potential.hpp"

namespace simol
{

  class DoubleWell : public Potential
  {
    public:
      DoubleWell(Input const& input);
      double operator()(double position) const;
      DVec gradient(double position) const;
      virtual double shiftToHarmonic() const;
      virtual DVec polynomialCoeffs() const;
    private:
      double height_;
      double interWell_;
      double center_;
  };

}

#endif