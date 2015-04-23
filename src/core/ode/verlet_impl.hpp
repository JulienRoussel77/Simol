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
    particle.speed_ -= timeStep * potential.derivative(particle.position_) / 2;
    particle.position_ += timeStep * particle.speed_ / particle.mass_;
    particle.speed_ -= timeStep * potential.derivative(particle.position_) / 2;
  }
}

#endif
