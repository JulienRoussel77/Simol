#ifndef POTENTIAL_IMPL_HPP
#define POTENTIAL_IMPL_HPP

#include <cmath>

namespace simol
{

  template<class ScalarType> inline
  Potential<ScalarType>::Potential(ScalarType parameter, ScalarType energy)
  :parameter_(parameter), energy_(energy)
  {}

  template<class ScalarType> inline
  ScalarType Potential<ScalarType>::operator()(ScalarType const & position) const
  { return parameter_*(1-std::cos(energy_*position)); }

  template<class ScalarType> inline
  ScalarType Potential<ScalarType>::derivative(ScalarType const & position) const
  { return parameter_*energy_*std::sin(energy_*position); }

  template<class ScalarType> inline
  ScalarType Potential<ScalarType>::force(ScalarType const & position) const
  { return -derivative(position); }
}


#endif
