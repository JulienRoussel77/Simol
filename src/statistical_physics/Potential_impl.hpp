#ifndef POTENTIAL_IMPL_HPP
#define POTENTIAL_IMPL_HPP

#include <cmath>

template<class ScalarType> inline
Potential<ScalarType>::Potential(ScalarType parameter, ScalarType pulsatance)
:parameter_(parameter), pulsatance_(pulsatance)
{}

template<class ScalarType> inline
ScalarType Potential<ScalarType>::operator()(ScalarType const & position) const
{ return parameter_*(1-std::cos(pulsatance_*position)); }

template<class ScalarType> inline
ScalarType Potential<ScalarType>::derivative(ScalarType const & position) const
{ return parameter_*pulsatance_*std::sin(pulsatance_*position); }




#endif
