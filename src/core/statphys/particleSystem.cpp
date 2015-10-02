#ifndef SIMOL_PARTICLESYSTEM_IMPL_HPP
#define SIMOL_PARTICLESYSTEM_IMPL_HPP

#include "particleSystem.hpp"

namespace simol
{
  
  ParticleSystem::ParticleSystem(Input const& input):
  currentTimeIteration_(0), 
  configuration_(input.numberOfParticles()),
  output(input.outputFilename())
  {
    for (size_t i = 0; i<input.numberOfParticles(); i++)
      configuration_[i] = Particle(input.mass(), input.initialPosition(i), input.initialSpeed(i));
  }
      
  Particle & ParticleSystem::particle(size_t index) 
  { return configuration_[index]; }

  std::vector< Particle > & ParticleSystem::configuration() 
  { return configuration_; }

  void ParticleSystem::simulate(double const timeStep, Dynamics * model)
  {
    //std::cout << "simulate !" << std::endl;
    for (auto&& particle : configuration_)
    {
      //verlet_scheme(particle,model->potential(),timeStep);

      model->update(particle, timeStep);
      output.display(currentTimeIteration_*timeStep, particle);
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
