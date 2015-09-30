#ifndef SIMOL_PARTICLESYSTEM_HPP
#define SIMOL_PARTICLESYSTEM_HPP

#include <vector>
#include "particle.hpp"
#include "dynamics.hpp"

namespace simol
{

  class ParticleSystem
  {
    typedef Particle ParticleType;

    public:

      ParticleSystem(size_t const numberOfParticles, 
                     double const & mass, 
                     double const & initialPosition, 
                     double const & initialSpeed);

      ParticleSystem(size_t const & numberOfParticles, 
                     double const & distance);

      ParticleType & particle(size_t index);

      std::vector<ParticleType> & configuration(); 

      void simulate(double const timeStep, 
                    Dynamics const* model, 
                    std::ofstream & outputFile);
      
      size_t size() const;
      
      

    private:
      
      size_t currentTimeIteration_;
      std::vector<ParticleType> configuration_;
  };

  // redémarrage : nombres aleatoires, juste pannes, tous les N pas de temps
  //

}






#include "particleSystem.ipp"

#endif
