#ifndef SIMOL_DOUBLEWELL_HPP
#define SIMOL_DOUBLEWELL_HPP

#include "Potential.hpp"

namespace simol
{

  class DoubleWell : public Potential
  {
    public:
      DoubleWell(Input const& input);
      virtual string classname() const {return "DoubleWell";}
      virtual double operator()(double position) const;
      virtual double symmetricValue(double position) const;
      virtual double skewsymmetricValue(double position) const;
      virtual double scalarGradient(double position) const;
      virtual double laplacian(double position) const;
      virtual double shiftToHarmonic() const;
      virtual DVec polyCoeffs() const;
    private:
      double interWell_;
      double coeff_;
  };

}

#endif