#ifndef SIMOL_PARTICLESYSTEM_IMPL_HPP
#define SIMOL_PARTICLESYSTEM_IMPL_HPP

#include "particleSystem.hpp"

namespace simol
{
  
  ParticleSystem* createSystem(Input  const& input, int const& indexOfReplica)
  {
    if (input.systemName() == "Isolated")
      return new Isolated(input, indexOfReplica);
    else if (input.systemName() == "BiChain")
      return new BiChain(input, indexOfReplica);
    else if (input.systemName() == "TriChain")
      return new TriChain(input, indexOfReplica);
    else
      std::cout << input.systemName() << " is not a valid system !" << std::endl;
    
    return 0;
  }
  
  ParticleSystem::ParticleSystem(Input const& input, int const& /*indexOfReplica*/):
  dimension_(input.dimension()),
  configuration_(input.numberOfParticles())
  {}
  
  Particle & ParticleSystem::getParticle(size_t index) 
  { return configuration_[index]; }
  
  const size_t& ParticleSystem::dimension() const
  {
		return dimension_;
	}

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
    for (size_t iOfIteration  =0; iOfIteration < model->numberOfIterations(); ++iOfIteration)
    {
      if ((10*iOfIteration) % model->numberOfIterations() == 0)
       cout << "---- Run " << (100 * iOfIteration) / model->numberOfIterations() << " % completed ----" << endl;
     
      
      //double instant = iOfIteration * model->timeStep();
      computeOutput(output, model, iOfIteration);
      writeOutput(output, iOfIteration);
      simulate(model);
    }
    computeFinalOutput(output, model);
    writeFinalOutput(output, model);
  }
  
  void ParticleSystem::simulate(Dynamics * model)
  {
    for (auto&& particle : configuration_)
      model->updateBefore(particle);
    
    computeAllForces(model);
    
    for (auto&& particle : configuration_)
      model->updateAfter(particle);
  }
  
  void ParticleSystem::computeOutput(Output& output, Dynamics const* model, size_t iOfIteration)
  {
    //if (output.verbose() > 0 && output.doOutput(iOfIteration))
    if (output.verbose() > 0)
    {
      output.kineticEnergy() = 0;
      output.potentialEnergy() = 0;
      //Calcul de la température et de l'énergie
      for (size_t iOfParticle = 0; iOfParticle < numberOfParticles(); iOfParticle++)
      {
				Particle& particle = configuration_[iOfParticle];
				//cout << "index = " << min(i+1, numberOfParticles()) << endl;
				//Particle& nextParticle = configuration_[min(i+1, numberOfParticles()-1)];
				output.kineticEnergy() += particle.kineticEnergy();
				output.potentialEnergy() += particle.potentialEnergy();
				//cout << i << " : " << particle.kineticEnergy() << " + " << particle.potentialEnergy() << endl;
				//cout << "momentum : " << particle.momentum() << endl;
			}
    }
    computeProfile(output, model, iOfIteration);
    model->updateAllControlVariates(output, configuration_, iOfIteration);
      
  }
  
  void ParticleSystem::writeOutput(Output& output, size_t iOfIteration)
  {
    if (output.verbose() > 0 && output.doOutput(iOfIteration) && iOfIteration >= 100)
      output.display(configuration_, iOfIteration);
  }
  
  void ParticleSystem::computeFinalOutput(Output& /*output*/, Dynamics const* /*model*/)
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
  
  
  
  
  
  
  BiChain::BiChain(Input const& input, int const& indexOfReplica):
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
  
    void BiChain::simulate(Dynamics * model)
  {
    //for (auto&& particle : configuration_)
      //model->updateBefore(particle);
    for (size_t i=0; i < numberOfParticles(); i++)
      model->updateBefore(getParticle(i));
    
    computeAllForces(model);
    
    assert(numberOfParticles() > 1);
    
    for (size_t i=0; i < numberOfParticles(); i++)
      model->updateAfter(getParticle(i));
    
    model->updateAfterLeft(getParticle(0));
    model->updateAfterRight(getParticle(numberOfParticles() - 1));
  }
  
  void BiChain::computeAllForces(Dynamics const* model)
  {
    //std::cout << "BiChain::computeAllForces" << std::endl;
    for (auto&& particle : configuration_)
      model->resetForce(particle);
    //for (auto&& particle : configuration_)
    //  model->computeForce(particle);
    model->interaction(ancorParticle_, configuration_[0]);
    for (size_t i = 0; i < numberOfParticles() - 1; i++)
      model->interaction(configuration_[i], configuration_[i+1]);
  }
  
  void BiChain::computeProfile(Output& output, Dynamics const* model, size_t iOfIteration)
  {
		for (size_t iOfParticle = 0; iOfParticle < numberOfParticles(); iOfParticle++)
    {
			Particle& particle = configuration_[iOfParticle];
			double flow = 0;
			if (iOfParticle == 0)
			{
				flow = model->gamma() * (model->temperatureLeft() - 2 * getParticle(0).kineticEnergy());
			}
			else
			{
				flow = - getParticle(iOfParticle).energyGrad(0) * getParticle(iOfParticle-1).momentum(0);
			}
			
			size_t midNumber = (numberOfParticles()-1) / 2;
			if (iOfParticle == midNumber)
				output.energyMidFlow() = flow;
			
			output.appendTemperatureProfile(2 * particle.kineticEnergy(), iOfIteration, iOfParticle);
			output.appendFlowProfile(flow, iOfIteration, iOfParticle);
		}
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
  
    void TriChain::simulate(Dynamics * model)
  {
    //for (auto&& particle : configuration_)
      //model->updateBefore(particle);
    for (size_t i=0; i < numberOfParticles(); i++)
      model->updateBefore(getParticle(i));
    
    computeAllForces(model);
    
    assert(numberOfParticles() > 1);
    
    for (size_t i=0; i < numberOfParticles(); i++)
      model->updateAfter(getParticle(i));
    
    model->updateAfterLeft(getParticle(0));
    model->updateAfterRight(getParticle(numberOfParticles() - 1));
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
    model->bending(configuration_[numberOfParticles() - 2], configuration_[numberOfParticles() - 1]);
  }
  
  void TriChain::computeProfile(Output& output, Dynamics const* model, size_t iOfIteration)
  {
		output.energySumFlow() = 0;
		for (size_t iOfParticle = 0; iOfParticle < numberOfParticles(); iOfParticle++)
    {
			//Particle& particle = configuration_[iOfParticle];
			double bending = 0;
			double flow = 0;
			
			if (iOfParticle == 0)
			{
				bending = ancorParticle2_.position(0) - 2 * getParticle(0).position(0) + getParticle(1).position(0);
				flow = model->gamma() * (model->temperatureLeft() - 2 * getParticle(0).kineticEnergy())
					- ancorParticle2_.energyGrad(0) * getParticle(0).momentum(0);
			}
			else if (iOfParticle < configuration_.size() - 1)
			{
				bending = getParticle(iOfParticle-1).position(0) - 2 * getParticle(iOfParticle).position(0) + configuration_[iOfParticle+1].position(0);
				flow = - getParticle(iOfParticle-1).energyGrad(0) * getParticle(iOfParticle).momentum(0) 
					+ getParticle(iOfParticle).energyGrad(0) * getParticle(iOfParticle-1).momentum(0);
				output.energySumFlow() += flow;
			}
			else
			{
				bending = 0;
				flow = - getParticle(iOfParticle-1).energyGrad(0) * getParticle(iOfParticle).momentum(0);
			}			
			
			size_t midNumber = (numberOfParticles()-1) / 2;
			if (iOfParticle == midNumber)
				output.energyMidFlow() = flow;
			
			output.appendBendingProfile(bending , iOfIteration, iOfParticle);
			output.appendTemperatureProfile(2 * getParticle(iOfParticle).kineticEnergy(), iOfIteration, iOfParticle);
			output.appendFlowProfile(flow, iOfIteration, iOfParticle);
		}
		output.energySumFlow() /= numberOfParticles()-2;
		//cout << output.energySumFlow() << endl;
	}
  
}

#endif
