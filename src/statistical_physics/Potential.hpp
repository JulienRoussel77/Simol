#ifndef POTENTIAL_HPP
#define POTENTIAL_HPP

template<class ScalarType>
class Potential
{
public:
  Potential(ScalarType parameter, ScalarType pulsatance);

  ScalarType operator()(ScalarType const & position) const;

  ScalarType derivative(ScalarType const & position) const;

private:
  ScalarType parameter_;
  ScalarType pulsatance_;
};

#include "Potential_impl.hpp"

#endif
