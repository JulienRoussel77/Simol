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
  configuration_(input.numberOfParticles())
  {}
  
  Particle & ParticleSystem::particle(size_t index) 
  { return configuration_[index]; }

  std::vector< Particle > & ParticleSystem::configuration() 
  { return configuration_; }
  
  size_t ParticleSystem::numberOfParticles() const
  {
    return configuration_.size();
  }
  
  void ParticleSystem::launch(Dynamics* model, Output& output)  
  {
    computeAllForces(model);
    for (size_t indexOfIteration  =0; indexOfIteration < model->numberOfIterations(); ++indexOfIteration)
    {
      if ((10*indexOfIteration) % model->numberOfIterations() == 0)
       cout << "---- Run " << (100 * indexOfIteration) / model->numberOfIterations() << " % completed ----" << endl;
     
      
      double instant = indexOfIteration * model->timeStep();
      computeOutput(output, model, indexOfIteration);
      writeOutput(output, indexOfIteration);
      simulate(model, output);
    }
    computeFinalOutput(output, model);
    writeFinalOutput(output, model);
  }
  
  void ParticleSystem::simulate(Dynamics * model, Output& output)
  {
    for (auto&& particle : configuration_)
      model->updateBefore(particle);
    
    computeAllForces(model);
    
    for (auto&& particle : configuration_)
      model->updateAfter(particle);
  }
  
  void ParticleSystem::computeOutput(Output& output, Dynamics const* model, size_t indexOfIteration)
  {
      
    //if (output.verbose() > 0 && output.doOutput(indexOfIteration))
    if (output.verbose() > 0)
    {
      output.kineticEnergy() = 0;
      output.potentialEnergy() = 0;
      //Calcul de la température et de l'énergie
      for (auto&& particle : configuration_)
      {
	output.kineticEnergy() += particle.kineticEnergy();
	output.potentialEnergy() += particle.potentialEnergy();
      }
      
      //Calcul de l'autocorrélation des vitesses
      /*if (output.doComputeCorrelations())
      {
	//output.appendAutocorrelationVelocity(particle(0).position(), indexOfIteration);
	output.appendAutocorrelationVelocity(particle(0).velocity(), indexOfIteration);
	output.appendAutocorrelationForce(particle(0).force() - model->externalForce(), indexOfIteration);
	
      }*/
     
    }
    model->updateAllControlVariates(output, configuration_, indexOfIteration);
      
  }
  
  void ParticleSystem::writeOutput(Output& output, size_t indexOfIteration)
  {
    if (output.verbose() > 0 && output.doOutput(indexOfIteration))
      output.display(configuration_, indexOfIteration);
  }
  
  void ParticleSystem::computeFinalOutput(Output& output, Dynamics const* model)
  {
    //output.responseForces() /= model->numberOfIterations();  
  }
  
  void ParticleSystem::writeFinalOutput(Output& output, Dynamics const* model)
  {
    //double time = model->timeStep() * model->numberOfIterations();
    
    output.finalDisplay(configuration_, model->externalForce(), model->numberOfIterations());
    if (output.doComputeCorrelations() &&  output.verbose() > 0)
      output.finalDisplayAutocorrelations();
  }
  
  Isolated::Isolated(Input const& input, int const& indexOfReplica):
  ParticleSystem(input, indexOfReplica)
  {
    for (size_t i = 0; i<input.numberOfParticles(); i++) 
      configuration_[i] = Particle(input.mass(), input.initialPosition(i), input.initialMomentum(i));
  }
      



  
  void Isolated::computeAllForces(Dynamics const* model)
  {
    //cout << "Isolated::computeAllForces" << endl;
    for (auto&& particle : configuration_)
    {
      model->resetForce(particle);
    }
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
  
    void Chain::simulate(Dynamics * model, Output& output)
  {
    //for (auto&& particle : configuration_)
      //model->updateBefore(particle);
    for (size_t i=1; i < numberOfParticles(); i++)
      model->updateBefore(particle(i));
    
    computeAllForces(model);
    
    assert(numberOfParticles() > 2);
    model->updateAfterLeft(particle(1));
    for (size_t i=1; i < numberOfParticles(); i++)
      model->updateAfter(particle(i));
    model->updateAfterRight(particle(numberOfParticles() - 1));
  }
  
  void Chain::computeAllForces(Dynamics const* model)
  {
    //std::cout << "Chain::computeAllForces" << std::endl;
    for (auto&& particle : configuration_)
      model->resetForce(particle);
    //for (auto&& particle : configuration_)
    //  model->computeForce(particle);
    for (size_t i = 0; i < numberOfParticles() - 1; i++)
      model->interaction(configuration_[i], configuration_[i+1]);
  }

}

#endif
