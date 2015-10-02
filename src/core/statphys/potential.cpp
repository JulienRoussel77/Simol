#ifndef POTENTIAL_IMPL_HPP
#define POTENTIAL_IMPL_HPP

#include "potential.hpp"

#include <cmath>

namespace simol
{

  Potential::Potential(double parameter, double energy)
  :parameter_(parameter), energy_(energy)
  {}
  
  Potential::Potential(Input const & input)
  :parameter_(input.potParameter()), energy_(2*M_PI/input.length())
  {}

  double Potential::operator()(dvec const & position) const
  { return parameter_*(1-std::cos(energy_*position(0))); }

  dvec Potential::derivative(dvec const & position) const
  { return parameter_*energy_*std::sin(energy_*position(0)); }

  dvec Potential::force(dvec const & position) const
  { return -derivative(position(0)); }
}


#endif
