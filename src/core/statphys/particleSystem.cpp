#ifndef SIMOL_PARTICLESYSTEM_IMPL_HPP
#define SIMOL_PARTICLESYSTEM_IMPL_HPP

#include "particleSystem.hpp"

namespace simol
{

  ParticleSystem::ParticleSystem(size_t const numberOfParticles, 
                                             std::vector<double> const & mass,
                                             std::vector<dvec> const & initialPosition,
                                             std::vector<dvec> const & initialSpeed)
  :currentTimeIteration_(0), 
   configuration_(numberOfParticles)
  {
    for (int i = 0; i<numberOfParticles; i++)
      configuration_[i] = Particle(mass[i], initialPosition[i], initialSpeed[i]);
  }
  
  ParticleSystem::ParticleSystem(Input const& input):
  currentTimeIteration_(0), 
  configuration_(input.numberOfParticles())
  {
    for (int i = 0; i<input.numberOfParticles(); i++)
      configuration_[i] = Particle(input.mass(), input.initialPosition(i), input.initialSpeed(i));
  }
      
  Particle & ParticleSystem::particle(size_t index) 
  { return configuration_[index]; }

  std::vector< Particle > & ParticleSystem::configuration() 
  { return configuration_; }

  void ParticleSystem::simulate(double const timeStep,
                                            Dynamics * model, 
                                            std::ofstream & outputFile)
  {
    for (auto&& particle : configuration_)
    {
      //verlet_scheme(particle,model->potential(),timeStep);
      model->update(particle, timeStep);
      //outputFile << particle.position();
      /*outputFile << currentTimeIteration_ * timeStep 
                 << " " << particle.position() 
                 << " " << particle.momentum() 
		 << " " << particle.kineticEnergy()
		 << " " << particle.potentialEnergy(model->potential())
		 << " " << particle.energy(model->potential())
                 << std::endl;*/
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
