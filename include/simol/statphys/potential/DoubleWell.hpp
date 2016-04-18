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
      Vector<double> gradient(double position) const;

    private:
      double height_;
      double interWell_;
  };

}

#endif