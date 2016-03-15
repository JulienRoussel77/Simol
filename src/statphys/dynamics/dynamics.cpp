#ifndef SIMOL_DYNAMICS_IPP
#define SIMOL_DYNAMICS_IPP

#include "dynamics.hpp"

using std::cout; 
using std::endl; 

namespace simol
{

	///
	///Calls the proper constructor of Dynamics and returns a pointer
  Dynamics* createDynamics(Input const& input, size_t iOfReplica)
  {
		//cout << "createDynamics(Input const& input, size_t iOfReplica)" << endl;
    if (input.dynamicsName() == "Hamiltonian")
      return new Hamiltonian(input, iOfReplica);
    else if (input.dynamicsName() == "Langevin")
      return new Langevin(input, iOfReplica);
    else if (input.dynamicsName() == "BoundaryLangevin")
      return new BoundaryLangevin(input, iOfReplica);
    else if (input.dynamicsName() == "Overdamped")
      return new Overdamped(input, iOfReplica);
    else
      std::cout << input.dynamicsName() << " is not a valid dynamics !" << std::endl;
    return 0;
  }
  
  ///
  ///Constructor
  Dynamics::Dynamics(Input const& input, int const& iOfReplica):
    timeStep_(input.timeStep(iOfReplica)), 
    nbOfIterations_(input.nbOfIterations(iOfReplica)), 
    nbOfThermalIterations_(input.nbOfThermalIterations(iOfReplica)), 
    nbOfBurningIterations_(input.nbOfBurningIterations(iOfReplica)), 
    externalForce_(input.dimension(), 0)
  {
    potential_ = createPotential(input);
    externalForce_(0) = input.externalForce(iOfReplica);
		
    cout << "externalForce = " << externalForce_(0) << endl;
    cout << "timeStep = " << timeStep() << endl;
    if (nbOfIterations_ < 1e6)
      cout << "nbOfIterations_ = " << nbOfIterations_ << endl;
    else
      cout << "nbOfIterations_ = " << nbOfIterations_/1e6 << "e6" << endl;
    cout << "finalTime = " << nbOfIterations_ * timeStep() << endl;
  }
  
  ///
  ///Destrucor
  Dynamics::~Dynamics()
  {
    delete potential_;
  }
  
  ///
  ///Set the Random Number Generator used by *this.
  void Dynamics::setRNG(RNG* rng) {rng_ = rng;}  
  ///
  ///Read-write accessor for the time step "dt"
  double& Dynamics::timeStep() {return timeStep_;}  
  ///
  ///Read-only accessor for the time step "dt"
  const double& Dynamics::timeStep() const {return timeStep_;}
  ///
  ///Read-write accessor for the number of iterations of the simulation
  size_t& Dynamics::nbOfIterations() {return nbOfIterations_;}
	///
  ///Read-only accessor for the number of iterations of the simulation
  const size_t& Dynamics::nbOfIterations() const {return nbOfIterations_;}
  ///
  ///Read-write accessor for the number of iterations of the thermalization
  size_t& Dynamics::nbOfThermalIterations() {return nbOfThermalIterations_;}
  ///
  ///Read-only accessor for the number of iterations of the thermalization
  const size_t& Dynamics::nbOfThermalIterations() const {return nbOfThermalIterations_;}
  ///
  ///Read-write accessor for the number of iterations of the burning
  size_t& Dynamics::nbOfBurningIterations() {return nbOfBurningIterations_;}
  ///
  ///Read-only accessor for the number of iterations of the burning
  const size_t& Dynamics::nbOfBurningIterations() const {return nbOfBurningIterations_;}
  ///
  ///Returns the total time "t_f" of the simulation
  double Dynamics::finalTime() const {return timeStep_ * nbOfIterations_;}
  ///
  ///Returns by value the potential of the dynamics
  Potential& Dynamics::potential() {return *potential_;}
  ///
  ///Evaluate the potential for the vector "position"
  double Dynamics::potential(Vector<double> const& position) const {return (*potential_)(position);}
  ///
  ///Evaluate the potential for the scalar "position"
  double Dynamics::potential(const double& distance) const {return (*potential_)(distance);}
  ///
  ///Evaluate the force for the scalar "position" (potential and external terms)
  Vector<double> Dynamics::force(Vector<double> const& position) const
  {
    return potential_->force(position) + externalForce_; 
  }
  ///
  ///Evaluate the laplacian of the potential for the vector "position"
  double Dynamics::laplacian(Vector<double> const& position) const
  {
    return potential_->laplacian(position); 
  }
  ///
  ///Read-only accessor for the external force
  const Vector<double>& Dynamics::externalForce() const {return externalForce_;}
  ///
  ///Write-read accessor for the external force
  Vector<double>& Dynamics::externalForce(){return externalForce_;}
  ///
  ///Read-only accessor for the i-th component of the external force
  const double& Dynamics::externalForce(int const& i) const {return externalForce_(i);}
  ///
  ///Write-read accessor for the i-th component of the external force
  double& Dynamics::externalForce(int const& i) {return externalForce_(i);}
	///
	///Returns the pointer of the Galerkin instance of *this
  Galerkin* Dynamics::galerkin() {return galerkin_;}
  
  
  
  ///
  ///Initialialize the momenta of the particles
  void Dynamics::initializeMomenta(vector<Particle>& /*configuration*/) {}
  ///
  ///Draw a momentum under the invariant measure at inverse temperature "localBeta"
  Vector<double> Dynamics::drawMomentum(double localBeta, double mass)
	{
		return sqrt(1 / (localBeta * mass)) * rng_->gaussian();
	}
	///
	///Draw a distance or a bending under the invariant measure at inverse temperature "localBeta"
  double Dynamics::drawPotLaw(double localBeta)
	{
		return potential_->drawLaw(localBeta, rng_);
	}
  ///Compute the mean distance or bending under the invariant measure
  ///Proceeds to a simple integral quadrature using rectangles
  double Dynamics::computeMeanPotLaw(double localBeta) const
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
  ///Set the force and the potential energy of a particle to zero
  void Dynamics::resetForce(Particle& particle) const
  {
    particle.potentialEnergy() = 0;
    particle.force().fill(0);
  }
	///
	///Computes the force and the energy of "particle" when the potential depends only on the positions
  void Dynamics::computeForce(Particle& particle) const
  {    
    particle.potentialEnergy() = potential(particle.position());
    particle.force() = force(particle.position());
  }
  
  ///Computes the force and the energy associated to this pair interaction, and updates these 2 fields
  ///The first 2 derivates of the potential are stored in "particle2"
  void Dynamics::interaction(Particle& particle1, Particle& particle2) const
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
  void Dynamics::triInteraction(Particle& particle1, Particle& particle2, Particle& particle3) const
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
  
  ///
  ///Standard first part of the numerical integration : half upadate on "p" and update on "q"
  void Dynamics::updateBefore(Particle& particle)
  {
    particle.momentum() += timeStep_ * particle.force() / 2;
    particle.position() += timeStep_ * particle.momentum() / particle.mass();
  }
  
  ///
  ///Standard second part of the numerical integration : half upadate on "p"
  void Dynamics::updateAfter(Particle& particle)
  {
    particle.momentum() += timeStep_ * particle.force() / 2;
  }
  
  ///
  ///Analytical integration of an Orstein-Uhlenbeck process of inverse T "localBeta"
  void Dynamics::updateOrsteinUhlenbeck(Particle& particle, double localBeta)
  {
    double alpha = exp(- gamma() / particle.mass() * timeStep_);    
    particle.momentum() = alpha * particle.momentum() + sqrt((1-pow(alpha, 2))/localBeta*particle.mass()) * rng_->gaussian();
  }
  
  ///
  ///Computes the quantities needed by the control variates (coefficients a, b, D) and {L \Phi}
  void Dynamics::updateAllControlVariates(Output& output, vector<Particle> const& configuration, size_t iOfIteration) const
  {
    Vector<double> q = configuration[0].position();
    Vector<double> p = configuration[0].momentum();
    Vector<double> qEnd = configuration[configuration.size()-1].position();
    VectorXd generatorOnBasis;
    generatorOnBasis = generatorOn(output.velocityCV(), configuration);
    output.velocityCV()->update(p(0), generatorOnBasis, configuration, iOfIteration);
    generatorOnBasis = generatorOn(output.forceCV(), configuration);
    output.forceCV()->update(potential_->derivative(q)(0), generatorOnBasis, configuration, iOfIteration);
    generatorOnBasis = generatorOn(output.lengthCV(), configuration);
    output.lengthCV()->update(qEnd(0), generatorOnBasis, configuration, iOfIteration);
    generatorOnBasis = generatorOn(output.midFlowCV(), configuration);
    output.midFlowCV()->update(output.energyMidFlow(), generatorOnBasis, configuration, iOfIteration);
		generatorOnBasis = generatorOn(output.sumFlowCV(), configuration);
    output.sumFlowCV()->update(output.energySumFlow(), generatorOnBasis, configuration, iOfIteration);
	}
  
  
  //#### Hamiltonian ####
  
  ///
  ///Constructs a Hamiltonian dynamics (constant energy)
  Hamiltonian::Hamiltonian(Input const& input, int const& iOfReplica):
		Dynamics(input, iOfReplica)
	{}
  
  ///
  ///Applie the generator of this dynamics to the basis functions of the CV
  MatrixXd Hamiltonian::generatorOn(ControlVariate const* controlVariate, vector<Particle> const& configuration) const
  {
    VectorXd result = VectorXd::Zero(controlVariate->nbOfFunctions());
    for (size_t iOfFunction=0; iOfFunction < controlVariate->nbOfFunctions(); iOfFunction++)
      for (size_t iOfParticle=0; iOfParticle < configuration.size(); iOfParticle++)
				result(iOfFunction) += configuration[iOfParticle].momentum().dot(controlVariate->gradientQ(configuration, iOfParticle, iOfFunction))
	      + force(configuration[iOfParticle].position()).dot(controlVariate->gradientP(configuration, iOfParticle, iOfFunction));  
    return result;      
  }
  
  
  
  //#### StochasticDynamics ####

  ///
  ///Constructs a purely virtual Dynamics for Dynamics using RNGs
  StochasticDynamics::StochasticDynamics(Input const& input, int const& iOfReplica):
    Dynamics(input, iOfReplica),
    xi_(input.xi())
  {
		cout << "xi = " << xi() << endl;
	}
  
  ///
  ///Read-only accessor for xi
  const double& StochasticDynamics::xi() const {return xi_;}
  ///
  ///Read-write accessor for xi
  double& StochasticDynamics::xi() {return xi_;}
	///
	///Returns the mean number of iterations between 2 random events
	int StochasticDynamics::xiNbOfIterations()
		{return 1 / (xi_ * timeStep_);}
	///
	///Returns true if the dynamics involves a Poisson process (momenta exchange)
	bool StochasticDynamics::doMomentaExchange() const
		{return xi_ > 0;}
	///
	///If the momenta exchange is activated, the times of future events are drawn
	void StochasticDynamics::initializeCountdown(Particle& particle)
	{
		if (doMomentaExchange())
			particle.countdown() = rng_->scalarExponential() * xiNbOfIterations(); // / (xi() * timestep());
		else
			particle.countdown() = -1;
	}
	///
	///Exchanges the momenta of the 2 particles if the time has come
	void StochasticDynamics::updateMomentaExchange(Particle& particle1, Particle& particle2)
	{
		if (particle2.countdown() == 0)
			{
				Vector<double> temp = particle2.momentum();
				particle2.momentum() = particle1.momentum();
				particle1.momentum() = temp;
				particle2.countdown() = rng_->scalarExponential() * xiNbOfIterations(); // / (xi() * timestep());
			}
			else
				particle2.countdown()--;
	}
     
  //#### UniformStochasticDynamics ####

  ///
  ///Constructs a purely virtual Dynamics for Dynamics with a global temperature
  UniformStochasticDynamics::UniformStochasticDynamics(Input const& input, int const& iOfReplica):
    StochasticDynamics(input, iOfReplica), 
    beta_(input.beta(iOfReplica)),
    temperature_(1/beta_)
  {}
  
  ///
  ///Initialize all the momenta under the invariant measure for inverse temperature "beta_"
  void UniformStochasticDynamics::initializeMomenta(vector<Particle>& configuration)
  {
    for (auto&& particle : configuration)
      //particle.momentum() = sqrt(1 / (beta() * particle.mass())) * rng_->gaussian();
			particle.momentum() = drawMomentum(beta(), particle.mass());
  }
  ///
  ///Returns the temperature
  double UniformStochasticDynamics::temperature() const {return temperature_;}
  ///
  ///Read-only access for the temperature
  const double& UniformStochasticDynamics::temperatureLeft() const {return temperature_;}
  ///
  ///Read-only access for the temperature
  const double& UniformStochasticDynamics::temperatureRight() const {return temperature_;}
  ///
  ///Read-only access for the inverse temperature
  const double& UniformStochasticDynamics::beta() const {return beta_;}
  ///
  ///Read-only access for the inverse temperature
  const double& UniformStochasticDynamics::betaLeft() const {return beta_;}
	///
  ///Read-only access for the inverse temperature
  const double& UniformStochasticDynamics::betaRight() const {return beta_;}


  //#### Langevin ####

  ///
  ///Constructor for the Langevin Dynamics with thermostats everywhere
  Langevin::Langevin(Input const& input, int const& iOfReplica):
    UniformStochasticDynamics(input, iOfReplica), 
    gamma_(input.gamma())
  {
		galerkin_ = createLangevinGalerkin(input);
	} 

  ///
  ///Read-only accessor for the intensity of the Orstein-Uhlenbeck process
  const double& Langevin::gamma() const {return gamma_;}
  
  ///
  ///Returns the amplitude of the brownian motion
  double Langevin::sigma() const
  {
    return sqrt(2*gamma_ / beta_);
  }

	///After refers to the fact that this step comes after updating the forces
	///Proceeds to the second half of the Verlet scheme, then integrate the O-U analytically
  void Langevin::updateAfter(Particle& particle)
  {
    particle.momentum() += timeStep_ * particle.force() / 2;
    double alpha = exp(- gamma_ / particle.mass() * timeStep_);    
    particle.momentum() = alpha * particle.momentum() + sqrt((1-pow(alpha, 2))/beta_*particle.mass()) * rng_->gaussian();
  } 
  
  ///Applies the generator of this dynamics to the basis functions of the CV
  ///Evaluate at the current state of "conifguration"
  MatrixXd Langevin::generatorOn(ControlVariate const* controlVariate, vector<Particle> const& configuration) const
  {
    VectorXd result = VectorXd::Zero(controlVariate->nbOfFunctions());
    for (size_t iOfFunction=0; iOfFunction < controlVariate->nbOfFunctions(); iOfFunction++)
      for (size_t iOfParticle=0; iOfParticle < configuration.size(); iOfParticle++)
			{
				result(iOfFunction) += dot(configuration[iOfParticle].momentum(), controlVariate->gradientQ(configuration, iOfParticle, iOfFunction))
					+ dot(configuration[iOfParticle].force(), controlVariate->gradientP(configuration, iOfParticle, iOfFunction))
					+ gamma_ * (- dot(configuration[iOfParticle].momentum(), controlVariate->gradientP(configuration, iOfParticle, iOfFunction))
							+ controlVariate->laplacianP(configuration, iOfParticle, iOfFunction) / beta_ );
			}
    return result;
  }
  
  ///
  ///Computes the quantities needed by the control variates (coefficients a, b, D) and {L \Phi}
  void Langevin::updateAllControlVariates(Output& output, vector<Particle> const& configuration, size_t iOfIteration) const
  {
		//cout << "Langevin::updateAllControlVariates(Output& output, vector<Particle> const& configuration, size_t iOfIteration)" << endl;
    Vector<double> q = configuration[0].position();
    Vector<double> p = configuration[0].momentum();
    VectorXd generatorOnBasis;
    generatorOnBasis = generatorOn(output.velocityCV(), configuration);
    output.velocityCV()->update(p(0), generatorOnBasis, configuration, iOfIteration);
		if (output.doOutput(iOfIteration))
			output.displayGeneratorOnBasis(output.outVelocitiesGenerator_, configuration, output.velocityCV(), iOfIteration*timeStep());
    
		/*generatorOnBasis = generatorOn(output.forceCV(), configuration);
    output.forceCV()->update(potential_->derivative(q)(0), generatorOnBasis, configuration, iOfIteration);
    generatorOnBasis = generatorOn(output.lengthCV(), configuration);
    output.lengthCV()->update(q(0), generatorOnBasis, configuration, iOfIteration);*/
		//cout << "end updateAllControlVariates" << endl;
  }
  
    
    
  //#### Overdamped ####

  ///
  ///Constructor for the Overdamped Langevin Dynamics
  Overdamped::Overdamped(Input const& input, int const& iOfReplica):
    UniformStochasticDynamics(input, iOfReplica)
  {}
  
  ///Before refers to the fact that this step comes before the forces update
  ///Do nothing
  void Overdamped::updateBefore(Particle& /*particle*/) {}
  
  // /!\ a changer dans le cas N != 1
  ///After refers to the fact that this step comes after the forces update
  ///Proceeds to a Euler-Maruyama scheme (in comment the second order is available)
  void Overdamped::updateAfter(Particle& particle)
  {
    particle.position() += timeStep_ * particle.force() + sqrt(2*timeStep_/beta_) * rng_->gaussian();
  
    //assert(particle.force()(0) == force(particle.position())(0));
    /*Vector<double> randomTerm = sqrt(2*timeStep_/beta_) * rng_->gaussian();
    Vector<double> qtilde = particle.position() + .5 * timeStep_ * particle.force() + .5 * randomTerm;
    particle.position() += timeStep_ * force(qtilde) + randomTerm;*/
  }
  
  ///Applies the generator of this dynamics to the basis functions of the CV
  ///Evaluate at the current state of "conifguration"
  MatrixXd Overdamped::generatorOn(ControlVariate const* controlVariate, vector<Particle> const& configuration) const
  {
    VectorXd result = VectorXd::Zero(controlVariate->nbOfFunctions());
    for (size_t iOfFunction=0; iOfFunction < controlVariate->nbOfFunctions(); iOfFunction++)
			for (size_t iOfParticle=0; iOfParticle < configuration.size(); iOfParticle++)
				result(iOfFunction) += controlVariate->laplacianQ(configuration, iOfParticle, iOfFunction) / beta_
					+ dot(force(configuration[iOfParticle].position()), controlVariate->gradientQ(configuration, iOfParticle, iOfFunction));
    return result;
  }

  
  //#### BoundaryLangevin ####

  ///
  ///Constructor for Langevin dynamics on chains, where there is a thermostat at each end
  BoundaryLangevin::BoundaryLangevin(Input const& input, int const& iOfReplica):
    StochasticDynamics(input, iOfReplica),
    betaLeft_(input.betaLeft()),
    betaRight_(input.betaRight()),
    temperatureLeft_(1/betaLeft_),
    temperatureRight_(1/betaRight()),
    gamma_(input.gamma()),
    tauBending_(input.tauBending())
  {
    cout << "deltaTemperature = " << deltaTemperature() << endl;    
  } 
  
  ///
  ///Read-only access for the inverse temperature at the left end
  const double& BoundaryLangevin::betaLeft() const {return betaLeft_;}
   ///
  ///Read-only access for the inverse temperature at the right end
  const double& BoundaryLangevin::betaRight() const {return betaRight_;}
  ///
  ///Returns the mean of the temperatures at the 2 ends
  double BoundaryLangevin::temperature() const
  {
    return (temperatureLeft_ + temperatureRight_) / 2;
  }
  ///
  ///Read-only access for the temperature at the left end
  const double& BoundaryLangevin::temperatureLeft() const {return temperatureLeft_;}
  ///
  ///Read-only access for the temperature at the right end
  const double& BoundaryLangevin::temperatureRight() const {return temperatureRight_;}
  
  ///Returns the difference between the mean temperature and the one at the left end
  ///This is equal to {eta}
  double BoundaryLangevin::deltaTemperature() const
  {
    return (temperatureLeft_ - temperatureRight_)/2;
  }
  ///
  ///Read-only accessor of the intensity of the O-U process
  const double& BoundaryLangevin::gamma() const {return gamma_;}
  ///
  ///Returns the amplitude of the brownian motion at the left end
  double BoundaryLangevin::sigmaLeft() const
  {
    return sqrt(2 * gamma_ / betaLeft_);
  }
  ///
  ///Returns the amplitude of the brownian motion at the right end
  double BoundaryLangevin::sigmaRight() const
  {
    return sqrt(2 * gamma_ / betaRight_);
  }
  ///
  ///Returns the bending constrain that is added on the right end of the chain
  const double& BoundaryLangevin::tauBending() const {return tauBending_;}
  ///
  ///Integrates the bending constraint on the (last) particle pair
  void BoundaryLangevin::bending(Particle& particle1, Particle& particle2) const
  {
    particle1.force(0) -= tauBending();
    particle2.force(0) += tauBending(); 
  }
  
  ///
  ///Initializes all the momenta, assuming an affine temperature profile
  void BoundaryLangevin::initializeMomenta(vector<Particle>& configuration)
  {
    for (size_t iOfParticle = 0; iOfParticle < configuration.size(); iOfParticle++)
    {
      Particle& particle = configuration[iOfParticle];
      double tempi_ = temperatureLeft_ + (iOfParticle + .5) * (temperatureRight_ - temperatureLeft_) / configuration.size();
      particle.momentum() = sqrt(tempi_ / particle.mass()) * rng_->gaussian();
    }
  }
  
  ///Applies the generator of this dynamics to the basis functions of the CV
  ///Evaluate at the current state of "conifguration"
  MatrixXd BoundaryLangevin::generatorOn(ControlVariate const* controlVariate, vector<Particle> const& configuration) const
  {
    size_t nbOfParticles = configuration.size();
    VectorXd result = VectorXd::Zero(controlVariate->nbOfFunctions());
    for (size_t iOfFunction=0; iOfFunction < controlVariate->nbOfFunctions(); iOfFunction++)
    {
      for (size_t iOfParticle=0; iOfParticle < configuration.size(); iOfParticle++)
      result(iOfFunction) += dot(configuration[iOfParticle].momentum(), controlVariate->gradientQ(configuration, iOfParticle, iOfFunction))
	      + dot(configuration[iOfParticle].force(), controlVariate->gradientP(configuration, iOfParticle, iOfFunction));
	      //if(false)	      
      result(iOfFunction) += gamma_ * (- dot(configuration[0].momentum(), controlVariate->gradientP(configuration, 0, iOfFunction))
			+ controlVariate->laplacianP(configuration, 0, iOfFunction) / betaLeft_
			- dot(configuration[nbOfParticles-1].momentum(), controlVariate->gradientP(configuration, nbOfParticles-1, iOfFunction))
			+ controlVariate->laplacianP(configuration, nbOfParticles-1, iOfFunction) / betaRight_);
    }
    return result;   
  }
 
    
  
}

#endif
