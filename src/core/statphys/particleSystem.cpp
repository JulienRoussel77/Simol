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
  currentTimeIteration_(0), 
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
    output.initialize(particle(0).momentum());
    computeAllForces(model);
    for (size_t instantIndex  =0; instantIndex < model->numberOfIterations(); ++instantIndex)
    {
      double instant = instantIndex * model->timeStep();
      computeOutput(output, model, model->timeStep()*instantIndex);
      writeOutput(output, model->timeStep()*instantIndex);
      simulate(model, output);

      ++currentTimeIteration_;
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
    
    ++currentTimeIteration_;
  }
  
  void ParticleSystem::computeOutput(Output& output, Dynamics const* model, double time)
  {
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
      
      //Calcul de la reponse moyenne à un forçage
      dvec averageForce(dimension_, 0);
      for (auto&& particle : configuration_)
	averageForce += particle.force() - model->externalForce();
      //positions.push_back(configuration_(0).position()); 
      averageForce /= numberOfParticles();
      output.responseForces() += averageForce;
      
      //Calcul de l'autocorrélation des vitesses
      if (output.doComputeCorrelations())
      {
	if (((int)(time/model->timeStep()) % output.decorrelationNumberOfIterations()) == 0)      
	{
	  output.timeRef() = time;
	  output.velocityRef() = particle(0).velocity();
	  output.forceRef() = particle(0).force() - model->externalForce();
	}
	//output.autocorrelationV.append(/model->timeStep(), output.velocityRef().dot(particle(0).velocity()));
	output.appendAutocorrelationV(particle(0).velocity(), time);
	output.appendAutocorrelationF(particle(0).force() - model->externalForce(), time);
	//cout << time/model->timeStep() << endl;
	//output.integratedAutocorrelationP() += model->timeStep() * output.refVelocity().dot(particle(0).momentum());
      }
    }
  }
  
  void ParticleSystem::writeOutput(Output& output, double time)
  {
    if (output.verbose() > 0)
      output.display(configuration_, time);
  }
  
  void ParticleSystem::computeFinalOutput(Output& output, Dynamics const* model)
  {
    //output.responseForces() /= model->numberOfIterations();  
  }
  
  void ParticleSystem::writeFinalOutput(Output& output, Dynamics const* model)
  {
    double time = model->timeStep() * model->numberOfIterations();
    
    output.finalDisplay(configuration_, model->externalForce(), model->numberOfIterations()*model->timeStep());
    if (output.doComputeCorrelations())
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
    for (auto&& particle : configuration_)
      model->updateBefore(particle);
    
    computeAllForces(model);
    
    assert(numberOfParticles() > 1);
    model->updateAfterLeft(particle(0));
    for (size_t i=1; i < numberOfParticles() - 1; i++)
      model->updateAfter(particle(i));
    model->updateAfterRight(particle(numberOfParticles() - 2));
    
    ++currentTimeIteration_;
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
