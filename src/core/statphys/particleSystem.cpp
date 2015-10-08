#ifndef SIMOL_PARTICLESYSTEM_IMPL_HPP
#define SIMOL_PARTICLESYSTEM_IMPL_HPP

#include "particleSystem.hpp"

using std::cout; 
using std::endl; 

namespace simol
{
  
  ParticleSystem* createSystem(Input  const& input, int const& indexOfReplica)
  {
    if (input.systemName() == "Isolated")
      return new Isolated(input, indexOfReplica);
    else if (input.systemName() == "Chain")
      return new Chain(input, indexOfReplica);
    else
      std::cout << input.systemName() << " is not a valid system !" << std::endl;
    
    return 0;
  }
  
  ParticleSystem::ParticleSystem(Input const& input, int const& indexOfReplica):
  dimension_(input.dimension()),
  numberOfParticles_(input.numberOfParticles()),
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
  
  void ParticleSystem::launch(Dynamics* model, Output& output, double const& timeStep, size_t const& numberOfIterations)  
  {
    output.initialize();
    computeAllForces(model);
    for (size_t instantIndex  =1; instantIndex < numberOfIterations; ++instantIndex)
    {
      double instant = instantIndex * timeStep;
      simulate(model, output, timeStep);
      computeOutput(output, model);
      writeOutput(output, timeStep*instantIndex);
      ++currentTimeIteration_;
    }
    computeFinalOutput(output, model, numberOfIterations);
    //writeFinalOutput(output, timeStep*numberOfIterations);
  }
  
  void ParticleSystem::computeOutput(Output& output, Dynamics* model)
  {
    dvec averageForce(dimension_, 0);
    for (auto&& particle : configuration_)
      averageForce += particle.force() - model->externalForce();
    //positions.push_back(configuration_(0).position()); 
    averageForce /= numberOfParticles_;
    output.sumForces() += averageForce;
  }
  
  void ParticleSystem::writeOutput(Output& output, double time)
  {
    if (output.verbose_ > 0)
      for (auto&& particle : configuration_)
	output.display(particle, time);
  }
  
  void ParticleSystem::computeFinalOutput(Output& output, Dynamics* model, size_t const& numberOfIterations)
  {
    output.responseForces().push_back(output.sumForces()/numberOfIterations);
  }
  
  void ParticleSystem::writeFinalOutput(Output& output, double time)
  {
    for (auto&& particle : configuration_)
      output.finalDisplay(particle, time);
  }
  
  Isolated::Isolated(Input const& input, int const& indexOfReplica):
  ParticleSystem(input, indexOfReplica)
  {
    for (size_t i = 0; i<input.numberOfParticles(); i++) 
      configuration_[i] = Particle(input.mass(), input.initialPosition(i), input.initialMomentum(i));
  }
      


  void Isolated::simulate(Dynamics * model, Output& output, double const& timeStep)
  {
    //std::cout << "simulate !" << std::endl;

    for (auto&& particle : configuration_)
      model->update(particle, timeStep);
    computeAllForces(model);
    //for (auto&& particle : configuration_)
    //  output.display(currentTimeIteration_*timeStep, particle);
  }
  
  void Isolated::computeAllForces(Dynamics const* model)
  {
    //cout << "Isolated::computeAllForces" << endl;
    for (auto&& particle : configuration_)
      model->resetForce(particle);
    for (auto&& particle : configuration_)
      model->computeForce(particle);
  }
  
  
  
    Chain::Chain(Input const& input, int const& indexOfReplica):
  ParticleSystem(input, indexOfReplica)
  {
    for (size_t i = 0; i<input.numberOfParticles(); i++) {
      configuration_[i] = Particle(input.mass(), input.initialPosition(i), input.initialMomentum(i));
      //std::cout << configuration_[i].force() << std::endl;
    }
  }
      


  void Chain::simulate(Dynamics* model, Output& output, double const& timeStep)
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
