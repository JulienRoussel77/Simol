#ifndef SIMOL_DYNAMICS_IPP
#define SIMOL_DYNAMICS_IPP

#include "dynamics.hpp"

using std::cout; 
using std::endl; 

namespace simol
{


  
  ///
  ///Constructor
  Dynamics::Dynamics(Input const& input):
    timeStep_(input.timeStep()), 
    nbOfIterations_(input.nbOfIterations()), 
    nbOfThermalIterations_(input.nbOfThermalIterations()), 
    nbOfBurningIterations_(input.nbOfBurningIterations()), 
    externalForce_(input.dimension(), 0),
    rng_(input.rng())
  {
    externalForce_(0) = input.externalForce();
		
    cout << "externalForce = " << externalForce_(0) << endl;
    cout << "timeStep = " << timeStep() << endl;
    if (nbOfIterations_ < 1e6)
      cout << "nbOfIterations_ = " << nbOfIterations_ << endl;
    else
      cout << "nbOfIterations_ = " << nbOfIterations_/1e6 << "e6" << endl;
    cout << "finalTime = " << nbOfIterations_ * timeStep() << endl;
  }
  

  
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
  ///Returns the total time "t_f" of the simulation
  double Dynamics::finalTime() const {return timeStep_ * nbOfIterations_;}
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
  
  const std::shared_ptr<RNG> Dynamics::rng() const {return rng_;}

  std::shared_ptr<RNG> Dynamics::rng() {return rng_;}
  
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
  ///Set the force and the potential energy of a particle to zero
  void Dynamics::resetForce(Particle& particle) const
  {
    particle.potentialEnergy() = 0;
    particle.force() = externalForce();
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
  

  
  
  //#### Hamiltonian ####
  
  ///
  ///Constructs a Hamiltonian dynamics (constant energy)
  Hamiltonian::Hamiltonian(Input const& input):
		Dynamics(input)
	{}
  

  
  
  
  //#### StochasticDynamics ####

  ///
  ///Constructs a purely virtual Dynamics for Dynamics using RNGs
  StochasticDynamics::StochasticDynamics(Input const& input):
    Dynamics(input),
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
  UniformStochasticDynamics::UniformStochasticDynamics(Input const& input):
    StochasticDynamics(input), 
    beta_(input.beta()),
    temperature_(1/beta_)
  {}
  

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
  Langevin::Langevin(Input const& input):
    UniformStochasticDynamics(input), 
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
  

  

  
    
    
  //#### Overdamped ####

  ///
  ///Constructor for the Overdamped Langevin Dynamics
  Overdamped::Overdamped(Input const& input):
    UniformStochasticDynamics(input)
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
  


  
  //#### BoundaryLangevin ####

  ///
  ///Constructor for Langevin dynamics on chains, where there is a thermostat at each end
  BoundaryLangevin::BoundaryLangevin(Input const& input):
    StochasticDynamics(input),
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
  

  

 
    
  
}

#endif
