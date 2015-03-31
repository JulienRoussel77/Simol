#ifndef SIMOL_VERLET_HPP
#define SIMOL_VERLET_HPP

#include "Particle_fwd.hpp"
#include "Potential_fwd.hpp"

namespace simol
{
  template<class ScalarType>
  void verlet(Particle<ScalarType> & particle, Potential<ScalarType> const & potential, double timeStep);
}

#include "verlet_impl.hpp"


#endif 
