#ifndef SIMOL_DYNAMICS_IPP
#define SIMOL_DYNAMICS_IPP

#include "dynamics.hpp"

using std::cout; 
using std::endl; 

namespace simol
{

    Dynamics* createDynamics(Input const& input, int const& indexOfReplica)
  {
    if (input.dynamicsName() == "Hamiltonian")
      return new Hamiltonian(input, indexOfReplica);
    else if (input.dynamicsName() == "Langevin")
      return new Langevin(input, indexOfReplica);
    else if (input.dynamicsName() == "BoundaryLangevin")
      return new BoundaryLangevin(input, indexOfReplica);
    else if (input.dynamicsName() == "Overdamped")
      return new Overdamped(input, indexOfReplica);
    else
      std::cout << input.dynamicsName() << " is not a valid dynamics !" << std::endl;
    return 0;
  }
  
  Dynamics::Dynamics(Input const& input, int const& indexOfReplica):timeStep_(input.timeStep()), numberOfIterations_(input.numberOfIterations()), externalForce_(input.dimension())
  {
    potential_ = createPotential(input);
    
    externalForce_(0) = input.externalForce(indexOfReplica);
    //cout << "externalForce_(0) = " << externalForce_(0) << endl;
  }
  
  Dynamics::~Dynamics()
  {
    delete potential_;
  }
  
  double& Dynamics::timeStep()
  {
    return timeStep_;
  }
  
  const double& Dynamics::timeStep() const
  {
    return timeStep_;
  }
  
  size_t& Dynamics::numberOfIterations()
  {
    return numberOfIterations_;
  }
  
  const size_t& Dynamics::numberOfIterations() const
  {
    return numberOfIterations_;
  }
  
  double Dynamics::finalTime() const
  {
    return timeStep_ * numberOfIterations_;
  }
  
  Potential const & Dynamics::potential() const
  { return *potential_; }
  
  double Dynamics::potential(dvec const& position) const
  {
    return (*potential_)(position); 
  }
  
  double Dynamics::potential(double const& distance) const
  {
    return (*potential_)(distance); 
  }
  
  dvec Dynamics::force(dvec const& position) const
  {
    return potential_->force(position) + externalForce_; 
  }
  
  const dvec& Dynamics::externalForce() const
  {
    return externalForce_; 
  }
  
    dvec& Dynamics::externalForce()
  {
    return externalForce_; 
  }
  
  const double& Dynamics::externalForce(int const& i) const
  {
    return externalForce_(i); 
  }
  
  double& Dynamics::externalForce(int const& i)
  {
    return externalForce_(i); 
  }
  
  void Dynamics::resetForce(Particle& particle) const
  {
    particle.potentialEnergy() = 0;
    particle.force().fill(0);
  }

  void Dynamics::computeForce(Particle& particle) const
  {    
    particle.potentialEnergy() = potential(particle.position());
    particle.force() = force(particle.position());
  }
  
  void Dynamics::interaction(Particle& particle1, Particle& particle2) const
  {
    dvec r12 = particle2.position() - particle1.position();
    //double d12 = r12.norm();
    double energy12 = potential(r12);
    dvec force12 = force(r12);
    
    particle1.potentialEnergy() += energy12 / 2;
    particle2.potentialEnergy() += energy12 / 2;
    particle1.force() += force12;
    particle2.force() -= force12;
  }
  
  void Dynamics::updateBefore(Particle& particle)
  {
    particle.momentum() += timeStep_ * particle.force() / 2;
    particle.position() += timeStep_ * particle.momentum() / particle.mass();
  }
  
  void Dynamics::updateAfter(Particle& particle)
  {
    particle.momentum() += timeStep_ * particle.force() / 2;
    particle.kineticEnergy() = pow(particle.momentum().norm(), 2) / particle.mass() / 2;
  }
  
  
  
  //#### Hamiltonian ####
  
  Hamiltonian::Hamiltonian(Input const& input, int const& indexOfReplica):Dynamics(input, indexOfReplica)
  {
  }
  
  
  
    //#### StochasticDynamics ####

  StochasticDynamics::StochasticDynamics(Input const& input, int const& indexOfReplica):
    Dynamics(input, indexOfReplica)
  {}
      
  void StochasticDynamics::setRNG(RNG* rng)
  {
    rng_ = rng;
  };
  
    
  //#### UniformStochasticDynamics ####

  UniformStochasticDynamics::UniformStochasticDynamics(Input const& input, int const& indexOfReplica):
    StochasticDynamics(input, indexOfReplica), 
    temperature_(input.temperature()), 
    beta_(1/temperature_)
  {}
  
  double const& UniformStochasticDynamics::temperature() const
  {
    return temperature_;
  }
    
  double const& UniformStochasticDynamics::beta() const
  {
    return beta_;
  }
  

  //#### Langevin ####

  Langevin::Langevin(Input const& input, int const& indexOfReplica):
    UniformStochasticDynamics(input, indexOfReplica), 
    gamma_(input.gamma())
  {} 

    
  double const& Langevin::gamma() const
  {
    return gamma_;
  }
  
  double Langevin::sigma() const
  {
    return sqrt(2*gamma_ / beta_);
  }


  void Langevin::updateAfter(Particle& particle)
  {
    particle.momentum() += timeStep_ * particle.force() / 2;
    double alpha = exp(- gamma_ / particle.mass() * timeStep_);    
    particle.momentum() = alpha * particle.momentum() + sqrt((1-pow(alpha, 2))/beta_*particle.mass()) * rng_->gaussian();

    particle.kineticEnergy() = pow(particle.momentum().norm(), 2) / particle.mass() / 2;
  } 
  
  
    
  //#### Overdamped ####

  Overdamped::Overdamped(Input const& input, int const& indexOfReplica):
    UniformStochasticDynamics(input, indexOfReplica)
  {}
  
    
  void Overdamped::updateBefore(Particle& particle)
  {
    particle.position() += timeStep_ * particle.force() + sqrt(2*timeStep_/beta_) * rng_->gaussian();
  }
  
  void Overdamped::updateAfter(Particle& particle){}
  
  
  
  //#### BoundaryLangevin ####

  BoundaryLangevin::BoundaryLangevin(Input const& input, int const& indexOfReplica):
    StochasticDynamics(input, indexOfReplica),
    betaLeft_(input.betaLeft()),
    betaRight_(input.betaRight()),
    temperatureLeft_(1/betaLeft_),
    temperatureRight_(1/betaRight()),
    gamma_(input.gamma())
  {} 
  
    double const& BoundaryLangevin::betaLeft() const
  {
    return betaLeft_;
  }
  
  double const& BoundaryLangevin::betaRight() const
  {
    return betaRight_;
  }
  
    double const& BoundaryLangevin::temperatureLeft() const
  {
    return temperatureLeft_;
  }
  
  double const& BoundaryLangevin::temperatureRight() const
  {
    return temperatureRight_;
  }

    
  double const& BoundaryLangevin::gamma() const
  {
    return gamma_;
  }
  
  double BoundaryLangevin::sigmaLeft() const
  {
    return sqrt(2 * gamma_ / betaLeft_);
  }
  
    double BoundaryLangevin::sigmaRight() const
  {
    return sqrt(2 * gamma_ / betaRight_);
  }
  
    void BoundaryLangevin::updateAfterLeft(Particle& particle)
  {
    particle.momentum() += timeStep_ * particle.force() / 2;
    double alpha = exp(- gamma_ / particle.mass() * timeStep_);    
    particle.momentum() = alpha * particle.momentum() + sqrt((1-pow(alpha, 2))/betaLeft_*particle.mass()) * rng_->gaussian();

    particle.kineticEnergy() = pow(particle.momentum().norm(), 2) / particle.mass() / 2;
  }
  
    void BoundaryLangevin::updateAfterRight(Particle& particle)
  {
    particle.momentum() += timeStep_ * particle.force() / 2;
    double alpha = exp(- gamma_ / particle.mass() * timeStep_);    
    particle.momentum() = alpha * particle.momentum() + sqrt((1-pow(alpha, 2))/betaRight_*particle.mass()) * rng_->gaussian();

    particle.kineticEnergy() = pow(particle.momentum().norm(), 2) / particle.mass() / 2;
  }
    
}

#endif
