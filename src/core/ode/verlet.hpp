#ifndef SIMOL_VERLET_HPP
#define SIMOL_VERLET_HPP

namespace simol 
{ 
  template<class ScalarType> 
  class Particle;
  
  template<class ScalarType> 
  class Potential;
}


namespace simol
{
  template<class ScalarType>
  void verlet(Particle<ScalarType> & particle, Potential<ScalarType> const & potential, double timeStep);
}

#include "verlet.ipp"


#endif 
