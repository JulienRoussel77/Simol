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
      { return configuration_[index]; }

      std::vector<ParticleType> & configuration() 
      { return configuration_; }

    private:
      ScalarType currentTime_;
      std::vector<ParticleType> configuration_;
  };

  // redémarrage : nombres aleatoires, juste pannes, tous les N pas de temps
  //

}

#include "ParticleSystem.ipp"

#endif
