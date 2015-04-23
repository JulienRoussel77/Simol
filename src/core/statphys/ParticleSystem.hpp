#ifndef SIMOL_PARTICLESYSTEM_HPP
#define SIMOL_PARTICLESYSTEM_HPP

#include <vector>
#include "Particle.hpp"

namespace simol
{

  template<class ScalarType>
  class ParticleSystem
  {
    typedef Particle<ScalarType> ParticleType;

    public:
      ParticleSystem(size_t const numberOfParticles, ScalarType const & mass, ScalarType const & initialPosition, ScalarType const & initialSpeed);

      ParticleType & particle(size_t index) 
      { return particles_[index]; }

      std::vector<ParticleType> & particles() 
      { return particles_; }

    private:
      ScalarType currentTime_;
      std::vector<ParticleType> particles_;
  };

  // redémarrage : nombres aleatoires, juste pannes, tous les N pas de temps
  //

}

#include "ParticleSystem_impl.hpp"

#endif
