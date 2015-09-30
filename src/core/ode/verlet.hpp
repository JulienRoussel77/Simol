#ifndef SIMOL_VERLET_HPP
#define SIMOL_VERLET_HPP

#include "particle.hpp"
#include "potential.hpp"

namespace simol
{
  //class Particle;
  class Potential;
  class HamiltonDynamics;
  
  void verlet(Particle & particle, 
              HamiltonDynamics const & model, 
              double timeStep)
  {
    particle.momentum_ += timeStep * model.potential().force(particle.position()) / 2;
    particle.position_ += timeStep * particle.momentum_ / particle.mass_;
    particle.momentum_ += timeStep * model.potential().force(particle.position()) / 2;
  }
  
  /*template<class double>
  void verlet(AtomChain const & model, 
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
