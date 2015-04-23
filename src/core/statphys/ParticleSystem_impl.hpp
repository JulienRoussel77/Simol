#ifndef SIMOL_PARTICLESYSTEM_IMPL_HPP
#define SIMOL_PARTICLESYSTEM_IMPL_HPP

namespace simol
{

  template<class ScalarType>
  ParticleSystem<ScalarType>::ParticleSystem(size_t const numberOfParticles, 
                                             ScalarType const & mass,
                                             ScalarType const & initialPosition,
                                             ScalarType const & initialSpeed)
  :currentTime_(0), 
   particles_(numberOfParticles, Particle<ScalarType>(mass,initialPosition,initialSpeed))
  {}

}

#endif
