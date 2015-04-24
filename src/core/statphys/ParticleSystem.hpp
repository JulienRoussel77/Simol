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

      void simulate(ScalarType const nextTime, Potential<ScalarType> const & potential, std::ofstream & outputFile)
      {
        ScalarType timeStep = nextTime - currentTime_;
        for (auto&& particle : configuration_)
        {
          verlet(particle,potential,timeStep);
          outputFile << currentTime_ 
                     << " " << particle.position() 
                     << " " << particle.speed() 
                     << std::endl;
        }
        currentTime_ = nextTime;

      }

    private:
      ScalarType currentTime_;
      std::vector<ParticleType> configuration_;
  };

  // redémarrage : nombres aleatoires, juste pannes, tous les N pas de temps
  //

}

#include "ParticleSystem.ipp"

#endif
