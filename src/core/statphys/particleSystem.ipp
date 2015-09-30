#ifndef SIMOL_PARTICLESYSTEM_IMPL_HPP
#define SIMOL_PARTICLESYSTEM_IMPL_HPP

namespace simol
{

  ParticleSystem::ParticleSystem(size_t const numberOfParticles, 
                                             double const & mass,
                                             double const & initialPosition,
                                             double const & initialSpeed)
  :currentTimeIteration_(0), 
   configuration_(numberOfParticles, Particle(mass,initialPosition,initialSpeed))
  {}
      
  Particle & ParticleSystem::particle(size_t index) 
  { return configuration_[index]; }

  std::vector< Particle > & ParticleSystem::configuration() 
  { return configuration_; }

  void ParticleSystem::simulate(double const timeStep,
                                            Dynamics const* model, 
                                            std::ofstream & outputFile)
  {
    for (auto&& particle : configuration_)
    {
      verlet_scheme(particle,model->potential(),timeStep);
      outputFile << currentTimeIteration_ * timeStep 
                 << " " << particle.position() 
                 << " " << particle.momentum() 
		 << " " << particle.kineticEnergy()
		 << " " << particle.potentialEnergy(model->potential())
		 << " " << particle.energy(model->potential())
                 << std::endl;
    }
    ++currentTimeIteration_;

  }
  
/*
  void ParticleSystem::simulate(double const timeStep,
                                            AtomChain const & model, 
                                            std::ofstream & outputFile)
  {
    for (auto & particle : configuration_)
      verlet(particle,model,timeStep);
    
    double alpha_first = std::exp(-gamma * timeStep / particle(0).mass());
    double alpha_last = std::exp(-gamma * timeStep / particle(N).mass())
    
    for(auto & particle : configuration_)
    {
      outputFile << currentTimeIteration_ * timeStep 
                 << " " << particle.position() 
                 << " " << particle.momentum() 
                 << std::endl;
    }
    ++currentTimeIteration_;

  }
*/

  size_t ParticleSystem::size() const
  {
    return configuration_.size();
  }
}

#endif
