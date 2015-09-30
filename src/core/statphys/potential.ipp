#ifndef POTENTIAL_IMPL_HPP
#define POTENTIAL_IMPL_HPP

#include <cmath>

namespace simol
{

  Potential::Potential(double parameter, double energy)
  :parameter_(parameter), energy_(energy)
  {}

  double Potential::operator()(double const & position) const
  { return parameter_*(1-std::cos(energy_*position)); }

  double Potential::derivative(double const & position) const
  { return parameter_*energy_*std::sin(energy_*position); }

  double Potential::force(double const & position) const
  { return -derivative(position); }
}


#endif
