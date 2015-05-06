#ifndef SIMOL_VERLET_HPP
#define SIMOL_VERLET_HPP

namespace simol 
{ 
  template<class ScalarType> 
  class Particle;
  
  template<class ScalarType> 
  class Potential;
  
  template<class ScalarType> 
  class HamiltonianDynamics;
}


namespace simol
{
  template<class ScalarType>
  void verlet(Particle<ScalarType> & particle, 
              HamiltonianDynamics<ScalarType> const & model, 
              double timeStep);
}

#include "verlet.ipp"


#endif 
