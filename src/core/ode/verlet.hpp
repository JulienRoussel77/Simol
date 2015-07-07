#ifndef SIMOL_VERLET_HPP
#define SIMOL_VERLET_HPP

#include "Particle.hpp"
#include "Potential.hpp"

namespace simol
{
  template<class ScalarType> class Particle;
  template<class ScalarType> class Potential;
  template<class ScalarType> class HamiltonDynamics;
  
  template<class ScalarType>
  void verlet(Particle<ScalarType> & particle, 
              HamiltonDynamics<ScalarType> const & model, 
              double timeStep)
  {
    particle.momentum_ += timeStep * model.potential().force(particle.position()) / 2;
    particle.position_ += timeStep * particle.momentum_ / particle.mass_;
    particle.momentum_ += timeStep * model.potential().force(particle.position()) / 2;
  }
  
  /*template<class ScalarType>
  void verlet(AtomChain<ScalarType> const & model, 
              double timeStep)
  {
    for( size_t index = 1; index < model.numberOfAtoms(); ++index )
    {
      model.momentum(index) += timeStep * model.potential().force(model.position(index)) / 2;
      model.position(index) += timeStep * model.momentum(index) / model.mass(index);
      model.momentum(index) += timeStep * model.potential().force(model.position(index)) / 2;
    }
  }*/

}



#endif 
