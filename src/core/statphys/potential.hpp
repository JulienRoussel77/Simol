#ifndef SIMOL_POTENTIAL_HPP
#define SIMOL_POTENTIAL_HPP

namespace simol
{

  class Potential
  {
    public:
      Potential(double parameter, double pulsatance);
      double operator()(double const & position) const;
      double derivative(double const & position) const;
      double force(double const & position) const;

    private:
      double parameter_;
      double energy_;
  };

}

#include "potential.ipp"

#endif
