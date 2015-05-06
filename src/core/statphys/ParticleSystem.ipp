#ifndef SIMOL_PARTICLESYSTEM_IMPL_HPP
#define SIMOL_PARTICLESYSTEM_IMPL_HPP

namespace simol
{

  template<class ScalarType>
  ParticleSystem<ScalarType>::ParticleSystem(size_t const numberOfParticles, 
                                             ScalarType const & mass,
                                             ScalarType const & initialPosition,
                                             ScalarType const & initialSpeed)
  :currentTimeIteration_(0), 
   configuration_(numberOfParticles, Particle<ScalarType>(mass,initialPosition,initialSpeed))
  {}
      
  template<class ScalarType>
  Particle<ScalarType> & ParticleSystem<ScalarType>::particle(size_t index) 
  { return configuration_[index]; }

  template<class ScalarType>
  std::vector< Particle<ScalarType> > & ParticleSystem<ScalarType>::configuration() 
  { return configuration_; }

  template<class ScalarType>
  void ParticleSystem<ScalarType>::simulate(ScalarType const timeStep,
                                            HamiltonDynamics<ScalarType> const & model, 
                                            std::ofstream & outputFile)
  {
    for (auto&& particle : configuration_)
    {
      verlet(particle,model,timeStep);
      outputFile << currentTimeIteration_ * timeStep 
                 << " " << particle.position() 
                 << " " << particle.momentum() 
                 << std::endl;
    }
    ++currentTimeIteration_;

  }

}

#endif
