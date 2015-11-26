#ifndef SIMOL_PARTICLESYSTEM_IMPL_HPP
#define SIMOL_PARTICLESYSTEM_IMPL_HPP

#include "particleSystem.hpp"

using std::cout;
using std::endl;
using std::min;
using std::max;

namespace simol
{
  
  ParticleSystem* createSystem(Input  const& input, int const& indexOfReplica)
  {
    if (input.systemName() == "Isolated")
      return new Isolated(input, indexOfReplica);
    else if (input.systemName() == "Chain")
      return new Chain(input, indexOfReplica);
    else if (input.systemName() == "TriChain")
      return new TriChain(input, indexOfReplica);
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
    model->initializeMomenta(configuration_);
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
      output.energyFlow() = 0;
      //Calcul de la température et de l'énergie
      for (size_t i = 0; i < numberOfParticles(); i++)
      {
	Particle& particle = configuration_[i];
	//cout << "index = " << min(i+1, numberOfParticles()) << endl;
	//Particle& nextParticle = configuration_[min(i+1, numberOfParticles()-1)];
	output.kineticEnergy() += particle.kineticEnergy();
	output.potentialEnergy() += particle.potentialEnergy();
	//cout << i << " : " << particle.kineticEnergy() << " + " << particle.potentialEnergy() << endl;
	//cout << "momentum : " << particle.momentum() << endl;
	//if (i != numberOfParticles() - 1)
	  //output.energyFlow() += - particle.energyGrad(0) * particle.momentum(0) / (numberOfParticles() - 1);
      }   
      //cout << output.kineticEnergy() + output.potentialEnergy() << endl;
      double inFlux = model->temperatureLeft() - 2*particle(0).kineticEnergy();
      double outFlux = 2*particle(numberOfParticles()-1).kineticEnergy()- model->temperatureRight();
      
      //output.energyFlow() = model->temperatureLeft() - 2*particle(0).kineticEnergy();
      //output.energyFlow() = 2*particle(numberOfParticles()-1).kineticEnergy() - model->temperatureRight();
      output.energyFlow() = (pow(model->temperatureRight(), 2) * inFlux + pow(model->temperatureLeft(), 2) * outFlux) / (pow(model->temperatureRight(), 2) + pow(model->temperatureLeft(), 2));
      //int midNumber = (numberOfParticles()-1) / 2;
      //output.energyFlow() = - particle(0).energyGrad(0) * particle(0).momentum(0);
      
      //output.energyFlow() = model->temperatureLeft() - 2*particle(0).kineticEnergy();

    }
    model->updateAllControlVariates(output, configuration_, indexOfIteration);
      
  }
  
  void ParticleSystem::writeOutput(Output& output, size_t indexOfIteration)
  {
    if (output.verbose() > 0 && output.doOutput(indexOfIteration) && indexOfIteration >= 100)
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
  ParticleSystem(input, indexOfReplica),
  ancorParticle_(input.dimension())
  {
    assert(configuration_.size() > 1);
    ancorParticle_.position(0) = 2 * input.initialPosition(0) - input.initialPosition(1);
    for (size_t i = 0; i<input.numberOfParticles(); i++) 
    {
      configuration_[i] = Particle(input.mass(), input.initialPosition(i), input.initialMomentum(i));
      //std::cout << configuration_[i].force() << std::endl;
    }
        configuration_[0].position(0) -= .1;

  }
  
    void Chain::simulate(Dynamics * model, Output& output)
  {
    //for (auto&& particle : configuration_)
      //model->updateBefore(particle);
    for (size_t i=0; i < numberOfParticles(); i++)
      model->updateBefore(particle(i));
    
    computeAllForces(model);
    
    assert(numberOfParticles() > 1);
    
    for (size_t i=0; i < numberOfParticles(); i++)
      model->updateAfter(particle(i));
    
    model->updateAfterLeft(particle(0));
    model->updateAfterRight(particle(numberOfParticles() - 1));
  }
  
  void Chain::computeAllForces(Dynamics const* model)
  {
    //std::cout << "Chain::computeAllForces" << std::endl;
    for (auto&& particle : configuration_)
      model->resetForce(particle);
    //for (auto&& particle : configuration_)
    //  model->computeForce(particle);
    model->interaction(ancorParticle_, configuration_[0]);
    for (size_t i = 0; i < numberOfParticles() - 1; i++)
      model->interaction(configuration_[i], configuration_[i+1]);
  }
  
  
  
  
  
  
  TriChain::TriChain(Input const& input, int const& indexOfReplica):
  ParticleSystem(input, indexOfReplica),
  ancorParticle1_(input.dimension()),
  ancorParticle2_(input.dimension())
  {
    assert(configuration_.size() > 1);
    ancorParticle1_.position(0) = 3 * input.initialPosition(0) - 2*input.initialPosition(1);
    ancorParticle2_.position(0) = 2 * input.initialPosition(0) - input.initialPosition(1);
    for (size_t i = 0; i<input.numberOfParticles(); i++) 
    {
      configuration_[i] = Particle(input.mass(), input.initialPosition(i), input.initialMomentum(i));
      //std::cout << configuration_[i].force() << std::endl;
    }
    configuration_[0].position(0) -= .1;
  }
  
    void TriChain::simulate(Dynamics * model, Output& output)
  {
    //for (auto&& particle : configuration_)
      //model->updateBefore(particle);
    for (size_t i=0; i < numberOfParticles(); i++)
      model->updateBefore(particle(i));
    
    computeAllForces(model);
    
    assert(numberOfParticles() > 1);
    
    for (size_t i=0; i < numberOfParticles(); i++)
      model->updateAfter(particle(i));
    
    model->updateAfterLeft(particle(0));
    model->updateAfterRight(particle(numberOfParticles() - 1));
  }
  
  void TriChain::computeAllForces(Dynamics const* model)
  {
    //std::cout << "TriChain::computeAllForces" << std::endl;
    for (auto&& particle : configuration_)
      model->resetForce(particle);
    //for (auto&& particle : configuration_)
    //  model->computeForce(particle);
    model->triInteraction(ancorParticle1_, ancorParticle2_, configuration_[0]);
    model->triInteraction(ancorParticle2_, configuration_[0], configuration_[1]);
    for (size_t i = 0; i < numberOfParticles() - 2; i++)
      model->triInteraction(configuration_[i], configuration_[i+1], configuration_[i+2]);
  }

}

#endif
