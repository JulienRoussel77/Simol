#ifndef SIMOL_PARTICLESYSTEM_IMPL_HPP
#define SIMOL_PARTICLESYSTEM_IMPL_HPP

#include "particleSystem.hpp"

namespace simol
{
  
  ParticleSystem* createSystem(Input  const& input)
  {
		cout << "createSystem(Input  const& input)" << endl;
    if (input.systemName() == "Isolated")
      return new Isolated(input);
    else if (input.systemName() == "BiChain")
      return new BiChain(input);
    else if (input.systemName() == "TriChain")
      return new TriChain(input);
    else
      std::cout << input.systemName() << " is not a valid system !" << std::endl;
    
    return 0;
  }
  
  ParticleSystem::ParticleSystem(Input const& input):
  dimension_(input.dimension()),
  configuration_(input.nbOfParticles()),
  settingsPath_(input.settingsPath()),
  rng_(input.rng())
  {
    potential_ = createPotential(input);
  }
  
  ///
  ///Destrucor
  ParticleSystem::~ParticleSystem()
  {
    delete potential_;
  }
  
  const Particle& ParticleSystem::getParticle(size_t index) const
  { return configuration_[index]; }
  
  Particle& ParticleSystem::getParticle(size_t index) 
  { return configuration_[index]; }
  
  const size_t& ParticleSystem::dimension() const
  {
		return dimension_;
	}
	
	const std::vector< Particle > & ParticleSystem::configuration() const
  { return configuration_; }

  std::vector< Particle > & ParticleSystem::configuration() 
  { return configuration_; }
  
  size_t ParticleSystem::nbOfParticles() const
  {
    return configuration_.size();
  }
  
  ///
  ///Returns by value the potential of the dynamics
  Potential& ParticleSystem::potential() {return *potential_;}
  ///
  ///Evaluate the potential for the vector "position"
  double ParticleSystem::potential(Vector<double> const& position) const {return (*potential_)(position);}
  ///
  ///Evaluate the potential for the scalar "position"
  double ParticleSystem::potential(const double& distance) const {return (*potential_)(distance);}
  ///
  ///Evaluate the force for the scalar "position" (potential and external terms)
  Vector<double> ParticleSystem::force(Vector<double> const& position) const
  {
    return potential_->force(position); 
  }
  ///
  ///Evaluate the laplacian of the potential for the vector "position"
  double ParticleSystem::laplacian(Vector<double> const& position) const
  {
    return potential_->laplacian(position); 
  }
  
    
  const std::shared_ptr<RNG> ParticleSystem::rng() const {return rng_;}

  std::shared_ptr<RNG> ParticleSystem::rng() {return rng_;}
   
  
  double ParticleSystem::boundaryPotEnergy() const
  {return 0;}
  
  ///
  ///Draw a momentum under the invariant measure at inverse temperature "localBeta"
  Vector<double> ParticleSystem::drawMomentum(double localBeta, double mass)
  {
    return sqrt(1 / (localBeta * mass)) * rng_->gaussian();
  }
  ///
  ///Draw a distance or a bending under the invariant measure at inverse temperature "localBeta"
  double ParticleSystem::drawPotLaw(double localBeta)
  {
    return potential_->drawLaw(localBeta, rng_);
  }
  ///Compute the mean distance or bending under the invariant measure
  ///Proceeds to a simple integral quadrature using rectangles
  double ParticleSystem::computeMeanPotLaw(double localBeta) const
  {
    double repFunc = 0;
    double qInteg = 0;
    size_t nbIntegrationNodes = 1000;
    double step = 8. / nbIntegrationNodes;
    Vector<double> deltaQ(1);
    for (size_t iOfNode = 0; iOfNode < nbIntegrationNodes; iOfNode++)
    {
      deltaQ(0) = - 4 + iOfNode * step;
      repFunc += exp(-localBeta * potential(deltaQ));
      qInteg += deltaQ(0) * exp(-localBeta * potential(deltaQ));
    }
    return qInteg / repFunc;
  }
  
  
  ///
  ///Computes the force and the energy of "particle" when the potential depends only on the positions
  void ParticleSystem::computeForce(Particle& particle) const
  {    
    particle.potentialEnergy() += potential(particle.position());
    particle.force() += force(particle.position());
  }
  
  ///Computes the force and the energy associated to this pair interaction, and updates these 2 fields
  ///The first 2 derivates of the potential are stored in "particle2"
  void ParticleSystem::interaction(Particle& particle1, Particle& particle2) const
  {
    Vector<double> r12 = particle2.position() - particle1.position();
    double energy12 = potential(r12);
    Vector<double> force12 = force(r12);    // = - v'(q_2 - q_1)
    double lapla12 = laplacian(r12);  // v"(q_2 - q_1)
    
    particle2.potentialEnergy() = energy12;
    particle1.force() -= force12;
    particle2.force() += force12;
    particle2.energyGrad() = -force12;    // v'(q_2 - q_1)
    particle2.energyLapla() = lapla12;    // v"(q_2 - q_1)
  }
  
  ///Computes the force and the energy associated to this triplet interaction, and updates these 2 fields
  ///The first 2 derivates of the potential are stored in "particle2"
  void ParticleSystem::triInteraction(Particle& particle1, Particle& particle2, Particle& particle3) const
  {
    Vector<double> delta = particle3.position() - 2*particle2.position() + particle1.position();
    //double d12 = r12.norm();
    double energy123 = potential(delta);
    Vector<double> force123 = force(delta);    // = - v'(r_2)
    double lapla123 = laplacian(delta);
    
    particle2.potentialEnergy() = energy123;
    particle1.force() += force123;
    particle2.force() -= 2*force123;
    particle3.force() += force123;
    particle2.energyGrad() = -force123;    // - v'(r_2)
    particle2.energyLapla() = lapla123;    // v''(r_2)
  }
  
  
  

  

  

  
  void ParticleSystem::writeOutput(Output& output, size_t iOfIteration)
  {
		if (output.verbose() > 0 && output.doOutput(iOfIteration))// && iOfIteration >= 100)
      output.display(configuration_, iOfIteration);
  }
  
  void ParticleSystem::computeFinalOutput(Output& /*output*/, Dynamics const& /*dyna*/) {}
  
  void ParticleSystem::writeFinalOutput(Output& output, Dynamics const& dyna)
  {    
    output.finalDisplay(configuration_, dyna.externalForce());
    if (output.doComputeCorrelations() &&  output.verbose() > 0)
      output.finalDisplayAutocorrelations();
  }
  
  //### Isolated ###
  
  Isolated::Isolated(Input const& input):
  ParticleSystem(input)
  {
    //for (size_t i = 0; i<input.nbOfParticles(); i++) 
    getParticle(0) = Particle(input.mass(), input.initialPosition(), input.initialMomentum());
		cout << "Particle initialized !" << endl;
		cout << getParticle(0).mass() << "  " << getParticle(0).position() << "  " << getParticle(0).momentum() << endl;
	}
      
  
  void Isolated::computeAllForces(Dynamics const& dyna)
  {
    dyna.resetForce(getParticle(0));
    computeForce(getParticle(0));
  }
  
  void Isolated::writeFinalOutput(Output& output, Dynamics const& dyna)
  {
    //double time = dyna.timeStep() * dyna.nbOfIterations();
    
    output.finalDisplay(configuration_, dyna.externalForce());
    if (output.doComputeCorrelations() &&  output.verbose() > 0)
      output.finalDisplayAutocorrelations();
		output.displayFinalVelocity(dyna.temperature(), dyna.externalForce(0), output.velocityCV_->nbOfFourier(), output.velocityCV_->nbOfHermite());
  }
  
  //### Fluid ###
  
  Fluid::Fluid(Input const& input):
  ParticleSystem(input)
  {
		assert(configuration_.size() > 1);
	  for (size_t i = 0; i<input.nbOfParticles(); i++) 
    {
      configuration_[i] = Particle(input.mass(), input.initialPosition(i), input.initialMomentum(i));
      //std::cout << configuration_[i].force() << std::endl;
    }
  }
	
	  void Fluid::computeAllForces(Dynamics const& dyna)
  {
    //std::cout << "Fluid::computeAllForces" << std::endl;
    for (auto&& particle : configuration_)
      dyna.resetForce(particle);
    //for (auto&& particle : configuration_)
    //  dyna.computeForce(particle);
    for (size_t i = 0; i < nbOfParticles()-1; i++)
			for (size_t j = i+1; j < nbOfParticles(); j++)
				interaction(configuration_[i], configuration_[j]);
  }
  
  void Fluid::writeFinalOutput(Output& output, Dynamics const& dyna)
  {
    //double time = dyna.timeStep() * dyna.nbOfIterations();
    
    output.finalDisplay(configuration_, dyna.externalForce());
    if (output.doComputeCorrelations() &&  output.verbose() > 0)
      output.finalDisplayAutocorrelations();
		//output.displayFinalFlow(dyna.temperature(), dyna.deltaTemperature());
  }
  
  //### Chain ###
  
  Chain::Chain(Input const& input):
  ParticleSystem(input)
  {
		assert(configuration_.size() > 1);
	  for (size_t i = 0; i<input.nbOfParticles(); i++) 
    {
      configuration_[i] = Particle(input.mass(), input.initialPosition(i), input.initialMomentum(i));
      //std::cout << configuration_[i].force() << std::endl;
    }
  }
  
  
  
  void Chain::computeAllForces(Dynamics const& model) {throw std::invalid_argument("Function undefined");};
  
  void Chain::thermalize(Dynamics& dyna)
  {
    //for (auto&& particle : configuration_)
      //dyna.updateBefore(particle);
    for (size_t i=0; i < nbOfParticles(); i++)
      dyna.updateBefore(getParticle(i));
    
    computeAllForces(dyna);
    
    for (size_t i=0; i < nbOfParticles(); i++)
      dyna.updateAfter(getParticle(i));
    
		for (size_t i=0; i < nbOfParticles(); i++)
		{
			double localTemperature = dyna.temperatureLeft() + i * dyna.deltaTemperature() / nbOfParticles();
			dyna.updateOrsteinUhlenbeck(getParticle(0), 1 / localTemperature);
		}
  }
  

	
	void Chain::computeProfile(Output& output, Dynamics const& model, size_t iOfIteration) const {throw std::invalid_argument("Function undefined");} 
  
  void Chain::writeFinalOutput(Output& output, Dynamics const& model) {throw std::invalid_argument("Function undefined");}
  
  //###### BiChain ######
  
  BiChain::BiChain(Input const& input):
  Chain(input),
  ancorParticle_(input.dimension())
  {
    //ancorParticle_.position(0) = 2 * input.initialPosition(0) - input.initialPosition(1);
		ancorParticle_.position(0) = 0;
  }
  

	

  

  
  void BiChain::computeAllForces(Dynamics const& dyna)
  {
    //std::cout << "BiChain::computeAllForces" << std::endl;
    for (auto&& particle : configuration_)
      dyna.resetForce(particle);
    //for (auto&& particle : configuration_)
    //  dyna.computeForce(particle);
    interaction(ancorParticle_, configuration_[0]);
    for (size_t i = 0; i < nbOfParticles() - 1; i++)
		{
      interaction(configuration_[i], configuration_[i+1]);
		}
  }
	
		void BiChain::computeProfile(Output& output, Dynamics const& dyna, size_t iOfIteration) const
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
				flow = dyna.gamma() * (dyna.temperatureLeft() - 2 * getParticle(0).kineticEnergy());
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
	
	void BiChain::writeFinalOutput(Output& output, Dynamics const& dyna)
  {
    //double time = dyna.timeStep() * dyna.nbOfIterations();
    
    output.finalDisplay(configuration_, dyna.externalForce());
    if (output.doComputeCorrelations() &&  output.verbose() > 0)
      output.finalDisplayAutocorrelations();
		output.displayFinalFlow(dyna.temperature(), dyna.deltaTemperature());
  }
  
  
  
  //###### TriChain ######
  
  
  TriChain::TriChain(Input const& input):
  Chain(input),
  ancorParticle1_(input.dimension()),
  ancorParticle2_(input.dimension())
  {
    ancorParticle1_.position(0) = 0;//3 * input.initialPosition(0) - 2*input.initialPosition(1);
    ancorParticle2_.position(0) = 0;//2 * input.initialPosition(0) - input.initialPosition(1);

  }
  
  
  void TriChain::computeAllForces(Dynamics const& dyna)
  {
    //std::cout << "TriChain::computeAllForces" << std::endl;
    for (auto&& particle : configuration_)
      dyna.resetForce(particle);
    //for (auto&& particle : configuration_)
    //  dyna.computeForce(particle);
    triInteraction(ancorParticle2_, ancorParticle1_, configuration_[0]);
    triInteraction(ancorParticle1_, configuration_[0], configuration_[1]);
    for (size_t i = 0; i < nbOfParticles() - 2; i++)
      triInteraction(configuration_[i], configuration_[i+1], configuration_[i+2]);
    dyna.bending(configuration_[nbOfParticles() - 2], configuration_[nbOfParticles() - 1]);
  }
  
  double TriChain::boundaryPotEnergy() const
  {return ancorParticle1_.potentialEnergy();}
  
  void TriChain::computeProfile(Output& output, Dynamics const& dyna, size_t iOfIteration) const
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
				flow = dyna.gamma() * (dyna.temperatureLeft() - 2 * getParticle(0).kineticEnergy())
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
	
	
	void TriChain::writeFinalOutput(Output& output, Dynamics const& dyna)
  {
    //double time = dyna.timeStep() * dyna.nbOfIterations();
    
    output.finalDisplay(configuration_, dyna.externalForce());
    if (output.doComputeCorrelations() &&  output.verbose() > 0)
      output.finalDisplayAutocorrelations();
		output.displayFinalFlow(dyna.temperature(), dyna.deltaTemperature(), dyna.tauBending(), dyna.xi());
  }
  
}

#endif
