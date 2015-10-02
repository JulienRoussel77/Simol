#ifndef SIMOL_POTENTIAL_HPP
#define SIMOL_POTENTIAL_HPP

#include "Vector.hpp"	
#include "input.hpp"

namespace simol
{

  class Potential
  {
    public:
      Potential(double parameter, double pulsatance);
      Potential(Input const& input);
      double operator()(dvec const & position) const;
      dvec derivative(dvec const & position) const;
      dvec force(dvec const & position) const;

    private:
      double parameter_;
      double energy_;
  };

}

//#include "potential.ipp"

#endif
