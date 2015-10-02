#ifndef SIMOL_PARTICLESYSTEM_HPP
#define SIMOL_PARTICLESYSTEM_HPP

#include <vector>
#include "particle.hpp"
#include "dynamics.hpp"
#include "output.hpp"

namespace simol
{

  class ParticleSystem
  {
    typedef Particle ParticleType;

    public:

      ParticleSystem(Input const& input);

      ParticleType & particle(size_t index);

      std::vector<ParticleType> & configuration(); 

      void simulate(double const timeStep, Dynamics * model);
      
      size_t size() const;
      
      

    private:
      
      size_t currentTimeIteration_;
      std::vector<ParticleType> configuration_;
      Output output;
  };

  // redémarrage : nombres aleatoires, juste pannes, tous les N pas de temps
  //

}






//#include "particleSystem.ipp"

#endif
