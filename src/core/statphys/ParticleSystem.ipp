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
   configuration_(numberOfParticles, Particle<ScalarType>(mass,initialPosition,initialSpeed))
  {}
      
  template<class ScalarType>
  Particle<ScalarType> & ParticleSystem<ScalarType>::particle(size_t index) 
  { return configuration_[index]; }

  template<class ScalarType>
  std::vector< Particle<ScalarType> > & ParticleSystem<ScalarType>::configuration() 
  { return configuration_; }

  template<class ScalarType>
  void ParticleSystem<ScalarType>::simulate(ScalarType const nextTime, 
                                            Potential<ScalarType> const & potential, 
                                            std::ofstream & outputFile)
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

}

#endif
