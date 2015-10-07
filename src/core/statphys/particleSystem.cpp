#ifndef SIMOL_PARTICLESYSTEM_IMPL_HPP
#define SIMOL_PARTICLESYSTEM_IMPL_HPP

#include "particleSystem.hpp"

using std::cout; 
using std::endl; 

namespace simol
{
  
  ParticleSystem* createSystem(Input  const& input)
  {
    if (input.systemName() == "Isolated")
      return new Isolated(input);
    else if (input.systemName() == "Chain")
      return new Chain(input);
    else
      std::cout << input.systemName() << " is not a valid system !" << std::endl;
    
    return 0;
  }
  
  ParticleSystem::ParticleSystem(Input const& input):
  currentTimeIteration_(0), 
  configuration_(input.numberOfParticles())
  {}
  
  Particle & ParticleSystem::particle(size_t index) 
  { return configuration_[index]; }

  std::vector< Particle > & ParticleSystem::configuration() 
  { return configuration_; }
  
  size_t ParticleSystem::size() const
  {
    return configuration_.size();
  }
  
  void ParticleSystem::launch(Dynamics* model, Output& output, double const& timeStep, int const& numberOfIterations)  
  {
    cout << "launch !" << endl;
    computeAllForces(model);
    for (size_t instantIndex  =1; instantIndex < numberOfIterations; ++instantIndex)
    {
      double instant = instantIndex * timeStep;
      simulate(model, output, timeStep, numberOfIterations);
      //computeOutput();
      //writeOutput();
    }
    writeOutput(output);
  }
  
  void ParticleSystem::computeOutput()
  {
    //positions.push_back(configuration_(0).position()); 
  }
  
  void ParticleSystem::writeOutput(Output& output)
  {
    for (auto&& particle : configuration_)
      output.display(particle);
  }
  
  Isolated::Isolated(Input const& input):
  ParticleSystem(input)
  {
    for (size_t i = 0; i<input.numberOfParticles(); i++) 
      configuration_[i] = Particle(input.mass(), input.initialPosition(i), input.initialSpeed(i));
  }
      


  void Isolated::simulate(Dynamics * model, Output& output, double const& timeStep, int const& numberOfIterations)
  {
    //std::cout << "simulate !" << std::endl;

    for (auto&& particle : configuration_)
      model->update(particle, timeStep);
    computeAllForces(model);
    //for (auto&& particle : configuration_)
    //  output.display(currentTimeIteration_*timeStep, particle);

    ++currentTimeIteration_;
  }
  
  void Isolated::computeAllForces(Dynamics const* model)
  {
    for (auto&& particle : configuration_)
      model->resetForce(particle);
    for (auto&& particle : configuration_)
      model->computeForce(particle);
  }
  
  
  
    Chain::Chain(Input const& input):
  ParticleSystem(input)
  {
    for (size_t i = 0; i<input.numberOfParticles(); i++) {
      configuration_[i] = Particle(input.mass(), input.initialPosition(i), input.initialSpeed(i));
      //std::cout << configuration_[i].force() << std::endl;
    }
  }
      


  void Chain::simulate(Dynamics* model, Output& output, double const& timeStep, int const& numberOfIterations)
  {
    //std::cout << "simulate !" << std::endl;

    for (auto&& particle : configuration_)
      model->update(particle, timeStep);
    computeAllForces(model);

    ++currentTimeIteration_;
  }
  
  void Chain::computeAllForces(Dynamics const* model)
  {
    //std::cout << "Chain::computeAllForces" << std::endl;
    for (auto&& particle : configuration_)
      model->resetForce(particle);
    //for (auto&& particle : configuration_)
    //  model->computeForce(particle);
    for (size_t i=0; i<size()-1; i++)
      model->interaction(configuration_[i], configuration_[i+1]);
  }

}

#endif
