#ifndef SIMOL_VERLET_IMPL_HPP
#define SIMOL_VERLET_IMPL_HPP

#include "Particle.hpp"
#include "Potential.hpp"

namespace simol
{
  template<class ScalarType> inline
  void verlet(Particle<ScalarType> & particle, 
              Potential<ScalarType> const & potential, 
              double timeStep)
  {
    particle.momentum_ += timeStep * potential.force(particle.position()) / 2;
    particle.position_ += timeStep * particle.momentum_ / particle.mass_;
    particle.momentum_ += timeStep * potential.force(particle.position()) / 2;
  }
}

#endif
