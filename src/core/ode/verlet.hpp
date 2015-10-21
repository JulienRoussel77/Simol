#ifndef SIMOL_VERLET_HPP
#define SIMOL_VERLET_HPP

#include "particle.hpp"
#include "potential.hpp"

/*namespace simol
{
  //class Particle;
  class Potential;
  class HamiltonDynamics;
  
  void verlet(Particle & particle, 
              HamiltonDynamics const & model, 
              double timeStep)
  {
    particle.momentum_ += timeStep * particle.force_ / 2;
    particle.position_ += timeStep * particle.momentum_ / particle.mass_;
    particle.momentum_ += timeStep * particle.force_ / 2;
    std::cout << particle.position_ << " " << particle.momentum_ << " " << particle.force_ << std::endl;
  }


}*/



#endif 
