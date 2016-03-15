#ifndef SIMOL_PARTICLESYSTEM_IMPL_HPP
#define SIMOL_PARTICLESYSTEM_IMPL_HPP

#include "particleSystem.hpp"

namespace simol
{
  
  ParticleSystem* createSystem(Input  const& input, int const& iOfReplica)
  {
		cout << "createSystem(Input  const& input, int const& iOfReplica)" << endl;
    if (input.systemName() == "Isolated")
      return new Isolated(input, iOfReplica);
    else if (input.systemName() == "BiChain")
      return new BiChain(input, iOfReplica);
    else if (input.systemName() == "TriChain")
      return new TriChain(input, iOfReplica);
    else
      std::cout << input.systemName() << " is not a valid system !" << std::endl;
    
    return 0;
  }
  
  ParticleSystem::ParticleSystem(Input const& input, int const& /*iOfReplica*/):
  dimension_(input.dimension()),
  configuration_(input.nbOfParticles()),
  settingsPath_(input.settingsPath())
  {}
  
  Particle & ParticleSystem::getParticle(size_t index) 
  { return configuration_[index]; }
  
  const size_t& ParticleSystem::dimension() const
  {
		return dimension_;
	}

  std::vector< Particle > & ParticleSystem::configuration() 
  { return configuration_; }
  
  size_t ParticleSystem::nbOfParticles() const
  {
    return configuration_.size();
  }
  
  double ParticleSystem::boundaryPotEnergy() const
  {return 0;}
  
  void ParticleSystem::launch(Dynamics* model, Output& output)  
  {
		cout << "Estimated time : " << 3.5 * nbOfParticles()/1024. * model->nbOfIterations() / 1e6 << " hours" << endl;
		//if (settingsPath_ == "")
		//	model->initializeMomenta(configuration_);
		//else
		//	model->startFrom(settingsPath_);
		
		initializeSystem(model);
		
    computeAllForces(model);
    for (size_t iOfIteration  =0; iOfIteration < model->nbOfIterations(); ++iOfIteration)
    {
      if ((10*iOfIteration) % model->nbOfIterations() == 0)
       cout << "---- Run " << (100 * iOfIteration) / model->nbOfIterations() << " % completed ----" << endl;
     
      
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
    if (output.verbose() > 0)
    {
      output.kineticEnergy() = 0;
      output.potentialEnergy() = 0;
      //Calcul de la température et de l'énergie
      for (size_t iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
      {
				Particle& particle = configuration_[iOfParticle];
				output.kineticEnergy() += particle.kineticEnergy();
				output.potentialEnergy() += particle.potentialEnergy();
			}
    }
    // In the case of the trichain we add the potential of the wall interaction
    output.potentialEnergy() += boundaryPotEnergy();
    computeProfile(output, model, iOfIteration);
    model->updateAllControlVariates(output, configuration_, iOfIteration);
  }
  
  void ParticleSystem::writeOutput(Output& output, size_t iOfIteration)
  {
		if (output.verbose() > 0 && output.doOutput(iOfIteration))// && iOfIteration >= 100)
      output.display(configuration_, iOfIteration);
  }
  
  void ParticleSystem::computeFinalOutput(Output& /*output*/, Dynamics const* /*model*/) {}
  
  void ParticleSystem::writeFinalOutput(Output& output, Dynamics const* model)
  {    
    output.finalDisplay(configuration_, model->externalForce());
    if (output.doComputeCorrelations() &&  output.verbose() > 0)
      output.finalDisplayAutocorrelations();
  }
  
  //### Isolated ###
  
  Isolated::Isolated(Input const& input, int const& iOfReplica):
  ParticleSystem(input, iOfReplica)
  {
    //for (size_t i = 0; i<input.nbOfParticles(); i++) 
    getParticle(0) = Particle(input.mass(), input.initialPosition(), input.initialMomentum());
		cout << "Particle initialized !" << endl;
		cout << getParticle(0).mass() << "  " << getParticle(0).position() << "  " << getParticle(0).momentum() << endl;
	}
      
  
  void Isolated::computeAllForces(Dynamics const* model)
  {
    model->resetForce(getParticle(0));
    model->computeForce(getParticle(0));
  }
  
  void Isolated::writeFinalOutput(Output& output, Dynamics const* model)
  {
    //double time = model->timeStep() * model->nbOfIterations();
    
    output.finalDisplay(configuration_, model->externalForce());
    if (output.doComputeCorrelations() &&  output.verbose() > 0)
      output.finalDisplayAutocorrelations();
		output.displayFinalVelocity(model->temperature(), model->externalForce(0), output.velocityCV_->nbOfFourier(), output.velocityCV_->nbOfHermite());
  }
  
  //### Fluid ###
  
  Fluid::Fluid(Input const& input, int const& iOfReplica):
  ParticleSystem(input, iOfReplica)
  {
		assert(configuration_.size() > 1);
	  for (size_t i = 0; i<input.nbOfParticles(); i++) 
    {
      configuration_[i] = Particle(input.mass(), input.initialPosition(i), input.initialMomentum(i));
      //std::cout << configuration_[i].force() << std::endl;
    }
  }
  
  void Fluid::simulate(Dynamics * model)
  {
    //for (auto&& particle : configuration_)
      //model->updateBefore(particle);
    for (size_t iOfParticle=0; iOfParticle < nbOfParticles(); iOfParticle++)
      model->updateBefore(getParticle(iOfParticle));
    
    computeAllForces(model);
    
    for (size_t iOfParticle=0; iOfParticle < nbOfParticles(); iOfParticle++)
      model->updateAfter(getParticle(iOfParticle));
	}
	
	  void Fluid::computeAllForces(Dynamics const* model)
  {
    //std::cout << "Fluid::computeAllForces" << std::endl;
    for (auto&& particle : configuration_)
      model->resetForce(particle);
    //for (auto&& particle : configuration_)
    //  model->computeForce(particle);
    for (size_t i = 0; i < nbOfParticles()-1; i++)
			for (size_t j = i+1; j < nbOfParticles(); j++)
				model->interaction(configuration_[i], configuration_[j]);
  }
  
  void Fluid::writeFinalOutput(Output& output, Dynamics const* model)
  {
    //double time = model->timeStep() * model->nbOfIterations();
    
    output.finalDisplay(configuration_, model->externalForce());
    if (output.doComputeCorrelations() &&  output.verbose() > 0)
      output.finalDisplayAutocorrelations();
		//output.displayFinalFlow(model->temperature(), model->deltaTemperature());
  }
  
  //### Chain ###
  
  Chain::Chain(Input const& input, int const& iOfReplica):
  ParticleSystem(input, iOfReplica)
  {
		assert(configuration_.size() > 1);
	  for (size_t i = 0; i<input.nbOfParticles(); i++) 
    {
      configuration_[i] = Particle(input.mass(), input.initialPosition(i), input.initialMomentum(i));
      //std::cout << configuration_[i].force() << std::endl;
    }
  }
  
  void Chain::thermalize(Dynamics * model)
  {
    //for (auto&& particle : configuration_)
      //model->updateBefore(particle);
    for (size_t i=0; i < nbOfParticles(); i++)
      model->updateBefore(getParticle(i));
    
    computeAllForces(model);
    
    for (size_t i=0; i < nbOfParticles(); i++)
      model->updateAfter(getParticle(i));
    
		for (size_t i=0; i < nbOfParticles(); i++)
		{
			double localTemperature = model->temperatureLeft() + i * model->deltaTemperature() / nbOfParticles();
			model->updateOrsteinUhlenbeck(getParticle(0), 1 / localTemperature);
		}
  }
  
  void Chain::simulate(Dynamics * model)
  {
    //for (auto&& particle : configuration_)
      //model->updateBefore(particle);
    for (size_t iOfParticle=0; iOfParticle < nbOfParticles(); iOfParticle++)
      model->updateBefore(getParticle(iOfParticle));
    
    computeAllForces(model);
    
    for (size_t iOfParticle=0; iOfParticle < nbOfParticles(); iOfParticle++)
      model->updateAfter(getParticle(iOfParticle));
    
    model->updateOrsteinUhlenbeck(getParticle(0), model->betaLeft());
    model->updateOrsteinUhlenbeck(getParticle(nbOfParticles() - 1), model->betaRight());
		
		if (model->doMomentaExchange())
			for (size_t iOfParticle=0; iOfParticle < nbOfParticles()-1; iOfParticle++)
				model->updateMomentaExchange(configuration_[iOfParticle], configuration_[iOfParticle+1]);
	}
  
  
  
  
  BiChain::BiChain(Input const& input, int const& iOfReplica):
  Chain(input, iOfReplica),
  ancorParticle_(input.dimension())
  {
    //ancorParticle_.position(0) = 2 * input.initialPosition(0) - input.initialPosition(1);
		ancorParticle_.position(0) = 0;
  }
  
  void BiChain::initializeSystem(Dynamics* model)
	{
		cout << "Initialization of the system...";cout.flush();
		double alpha, localTemp, localDist;
		//ofstream outTest("test.txt");
		for (size_t iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
		{
			alpha = iOfParticle / (double) nbOfParticles();			
			localTemp = (1-alpha) * model->temperatureLeft() + alpha * model->temperatureRight();
			//localBending = model->computeMeanPotLaw(1/localTemp);
			localDist = model->drawPotLaw(1/localTemp);
			//outTest << iOfParticle << " " << localBending << endl;
			//cout << "bending = " << localBending << " / mean = " << model->computeMeanPotLaw(1/localTemp) << endl;
			
			double prevPosition = (iOfParticle>0)?getParticle(iOfParticle-1).position(0):0;
			
			
			getParticle(iOfParticle).position(0) = prevPosition + localDist;
			//cout << -position2 << " + 2 * " << position1 << " + " << localBending << " = " << getParticle(iOfParticle).position(0) << endl;
			getParticle(iOfParticle).momentum() = model->drawMomentum(1/localTemp, getParticle(iOfParticle).mass());
		
			model->initializeCountdown(getParticle(iOfParticle));
		}
		
		cout << "Done ! / Thermalization...";cout.flush();
		
		for (size_t iOfIteration  =0; iOfIteration < model->nbOfThermalIterations(); ++iOfIteration)
    {
			thermalize(model);
		}
		
		cout << "Done ! / Burning...";cout.flush();
		
		for (size_t iOfIteration  =0; iOfIteration < model->nbOfBurningIterations(); ++iOfIteration)
    {
			simulate(model);
		}
		cout << "Done !" << endl;
	}
	

  

  
  void BiChain::computeAllForces(Dynamics const* model)
  {
    //std::cout << "BiChain::computeAllForces" << std::endl;
    for (auto&& particle : configuration_)
      model->resetForce(particle);
    //for (auto&& particle : configuration_)
    //  model->computeForce(particle);
    model->interaction(ancorParticle_, configuration_[0]);
    for (size_t i = 0; i < nbOfParticles() - 1; i++)
		{
      model->interaction(configuration_[i], configuration_[i+1]);
		}
  }
	
		void BiChain::computeProfile(Output& output, Dynamics const* model, size_t iOfIteration)
  {
		output.energySumFlow() = 0;
		size_t midNb = (nbOfParticles()-1) / 2;
		for (size_t iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
    {
			//Particle& particle = configuration_[iOfParticle];
			double dist = 0;
			double flow = 0;
			double potTempTop = 0;
			double potTempBot = 0;
			
			if (iOfParticle == 0)
			{
				dist = getParticle(0).position(0) - ancorParticle_.position(0);
				flow = model->gamma() * (model->temperatureLeft() - 2 * getParticle(0).kineticEnergy());
			}
			else
			{
				// dist is r_iOfParticle
				dist = getParticle(iOfParticle).position(0) - getParticle(iOfParticle-1).position(0);
				// flow is j_iOfParticle
				flow = - getParticle(iOfParticle).energyGrad(0) * getParticle(iOfParticle-1).momentum(0);
				if (iOfParticle != nbOfParticles()-1)
				{
					output.energySumFlow() += flow;					
					if (iOfParticle == midNb)
						output.energyMidFlow() = flow;
				}
			}
			
			if (iOfParticle == nbOfParticles()-1)
			{
				potTempTop = pow(getParticle(iOfParticle).energyGrad(0), 2);
				potTempBot = getParticle(iOfParticle).energyLapla();
			}
			else
			{
				potTempTop = pow(getParticle(iOfParticle).energyGrad(0) - getParticle(iOfParticle+1).energyGrad(0), 2);
				potTempBot = getParticle(iOfParticle).energyLapla() + getParticle(iOfParticle+1).energyLapla();
			}
			
			output.appendBendistProfile(dist , iOfIteration, iOfParticle);
			output.appendKinTempProfile(2 * getParticle(iOfParticle).kineticEnergy(), iOfIteration, iOfParticle);
			output.appendPotTempTopProfile(potTempTop, iOfIteration, iOfParticle);
			output.appendPotTempBotProfile(potTempBot, iOfIteration, iOfParticle);
			output.appendFlowProfile(flow, iOfIteration, iOfParticle);
		}
		output.energySumFlow() /= (nbOfParticles()-2.);
		//cout << output.energySumFlow() << endl;
	}
	
	void BiChain::writeFinalOutput(Output& output, Dynamics const* model)
  {
    //double time = model->timeStep() * model->nbOfIterations();
    
    output.finalDisplay(configuration_, model->externalForce());
    if (output.doComputeCorrelations() &&  output.verbose() > 0)
      output.finalDisplayAutocorrelations();
		output.displayFinalFlow(model->temperature(), model->deltaTemperature());
  }
  
  
  
  
  
  
  TriChain::TriChain(Input const& input, int const& iOfReplica):
  Chain(input, iOfReplica),
  ancorParticle1_(input.dimension()),
  ancorParticle2_(input.dimension())
  {
    ancorParticle1_.position(0) = 0;//3 * input.initialPosition(0) - 2*input.initialPosition(1);
    ancorParticle2_.position(0) = 0;//2 * input.initialPosition(0) - input.initialPosition(1);

  }
  
  void TriChain::initializeSystem(Dynamics* model)
	{
		cout << "Initialization of the system...";cout.flush();
		double alpha, localTemp, localBending;
		//ofstream outTest("test.txt");
		for (size_t iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
		{
			alpha = iOfParticle / (double) nbOfParticles();			
			localTemp = (1-alpha) * model->temperatureLeft() + alpha * model->temperatureRight();
			//localBending = model->computeMeanPotLaw(1/localTemp);
			localBending = model->drawPotLaw(1/localTemp);
			//outTest << iOfParticle << " " << localBending << endl;
			//cout << "bending = " << localBending << " / mean = " << model->computeMeanPotLaw(1/localTemp) << endl;
			
			double position1 = (iOfParticle>0)?getParticle(iOfParticle-1).position(0):0;
			double position2 = (iOfParticle>1)?getParticle(iOfParticle-2).position(0):0;
			
			
			getParticle(iOfParticle).position(0) = -position2 + 2 * position1 + localBending;
			//cout << -position2 << " + 2 * " << position1 << " + " << localBending << " = " << getParticle(iOfParticle).position(0) << endl;
			getParticle(iOfParticle).momentum() = model->drawMomentum(1/localTemp, getParticle(iOfParticle).mass());
			
			model->initializeCountdown(getParticle(iOfParticle));
		}
		
		cout << "Done ! / Thermalization...";cout.flush();
		
		for (size_t iOfIteration  =0; iOfIteration < model->nbOfThermalIterations(); ++iOfIteration)
    {
			thermalize(model);
		}
		
		cout << "Done ! / Burning...";cout.flush();
		
		for (size_t iOfIteration  =0; iOfIteration < model->nbOfBurningIterations(); ++iOfIteration)
    {
			simulate(model);
		}
		cout << "Done !" << endl;
	}
  
  void TriChain::computeAllForces(Dynamics const* model)
  {
    //std::cout << "TriChain::computeAllForces" << std::endl;
    for (auto&& particle : configuration_)
      model->resetForce(particle);
    //for (auto&& particle : configuration_)
    //  model->computeForce(particle);
    model->triInteraction(ancorParticle2_, ancorParticle1_, configuration_[0]);
    model->triInteraction(ancorParticle1_, configuration_[0], configuration_[1]);
    for (size_t i = 0; i < nbOfParticles() - 2; i++)
      model->triInteraction(configuration_[i], configuration_[i+1], configuration_[i+2]);
    model->bending(configuration_[nbOfParticles() - 2], configuration_[nbOfParticles() - 1]);
  }
  
  double TriChain::boundaryPotEnergy() const
  {return ancorParticle1_.potentialEnergy();}
  
  void TriChain::computeProfile(Output& output, Dynamics const* model, size_t iOfIteration)
  {
		output.energySumFlow() = 0;
		for (size_t iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
    {
			//Particle& particle = configuration_[iOfParticle];
			double bending = 0;
			double flow = 0;
			double potTempTop = 0;
			double potTempBot = 0;
			
			if (iOfParticle == 0)
			{
				bending = ancorParticle2_.position(0) - 2 * getParticle(0).position(0) + getParticle(1).position(0);
				flow = model->gamma() * (model->temperatureLeft() - 2 * getParticle(0).kineticEnergy())
					- ancorParticle2_.energyGrad(0) * getParticle(0).momentum(0);
				potTempTop = pow(- ancorParticle2_.energyGrad(0) + 2*getParticle(0).energyGrad(0) - getParticle(1).energyGrad(0), 2);
				potTempBot = ancorParticle2_.energyLapla() + 4*getParticle(0).energyLapla() + getParticle(1).energyLapla();
				
			}
			else if (iOfParticle < configuration_.size() - 1)
			{
				// bending is k_iOfParticle
				bending = getParticle(iOfParticle-1).position(0) - 2 * getParticle(iOfParticle).position(0) + configuration_[iOfParticle+1].position(0);
				// flow is j_iOfParticle
				flow = - getParticle(iOfParticle-1).energyGrad(0) * getParticle(iOfParticle).momentum(0) 
					+ getParticle(iOfParticle).energyGrad(0) * getParticle(iOfParticle-1).momentum(0);

				potTempTop = pow(- getParticle(iOfParticle-1).energyGrad(0) + 2*getParticle(iOfParticle).energyGrad(0) - getParticle(iOfParticle+1).energyGrad(0), 2);
				potTempBot = getParticle(iOfParticle-1).energyLapla() + 4*getParticle(iOfParticle).energyLapla() + getParticle(iOfParticle+1).energyLapla();
				output.energySumFlow() += flow;
			}
			else
			{
				bending = nan("");
				flow = - getParticle(iOfParticle-1).energyGrad(0) * getParticle(iOfParticle).momentum(0);
				potTempTop = pow(- getParticle(iOfParticle-1).energyGrad(0), 2);
				potTempBot = getParticle(iOfParticle-1).energyLapla();
			}			
			
			size_t midNb = (nbOfParticles()-1) / 2;
			if (iOfParticle == midNb)
				output.energyMidFlow() = flow;
			
			output.appendBendistProfile(bending , iOfIteration, iOfParticle);
			output.appendKinTempProfile(2 * getParticle(iOfParticle).kineticEnergy(), iOfIteration, iOfParticle);
			output.appendPotTempTopProfile(potTempTop, iOfIteration, iOfParticle);
			output.appendPotTempBotProfile(potTempBot, iOfIteration, iOfParticle);
			output.appendFlowProfile(flow, iOfIteration, iOfParticle);
		}
		output.energySumFlow() /= (nbOfParticles()-2.);
		//cout << output.energySumFlow() << endl;
	}
	
	
	void TriChain::writeFinalOutput(Output& output, Dynamics const* model)
  {
    //double time = model->timeStep() * model->nbOfIterations();
    
    output.finalDisplay(configuration_, model->externalForce());
    if (output.doComputeCorrelations() &&  output.verbose() > 0)
      output.finalDisplayAutocorrelations();
		output.displayFinalFlow(model->temperature(), model->deltaTemperature(), model->tauBending(), model->xi());
  }
  
}

#endif
