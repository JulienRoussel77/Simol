#ifndef SIMOL_DYNAMICS_IPP
#define SIMOL_DYNAMICS_IPP

#include "dynamics.hpp"

using std::cout; 
using std::endl; 

namespace simol
{

    Dynamics* createDynamics(Input const& input, size_t indexOfReplica)
  {
    if (input.dynamicsName() == "Hamiltonian")
      return new Hamiltonian(input, indexOfReplica);
    else if (input.dynamicsName() == "Langevin")
      return new Langevin(input, indexOfReplica);
    else if (input.dynamicsName() == "BoundaryLangevin")
      return new BoundaryLangevin(input, indexOfReplica);
    else if (input.dynamicsName() == "Overdamped")
      return new Overdamped(input, indexOfReplica);
    else if (input.dynamicsName() == "BoundaryOverdamped")
      return new BoundaryOverdamped(input, indexOfReplica);
    else
      std::cout << input.dynamicsName() << " is not a valid dynamics !" << std::endl;
    return 0;
  }
  
  Dynamics::Dynamics(Input const& input, int const& indexOfReplica):
    timeStep_(input.timeStep(indexOfReplica)), 
    numberOfIterations_(input.numberOfIterations(indexOfReplica)), 
    externalForce_(input.dimension())
  {
    potential_ = createPotential(input);
    externalForce_(0) = input.externalForce(indexOfReplica);
    cout << "externalForce = " << externalForce_(0) << endl;
    cout << "timeStep = " << timeStep() << endl;
    cout << "numberOfIterations_ = " << numberOfIterations_ << endl;
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
  
  Potential& Dynamics::potential()
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
    //cout << "Dynamics::updateBefore" <<endl;
    particle.momentum() += timeStep_ * particle.force() / 2;
    particle.position() += timeStep_ * particle.momentum() / particle.mass();
  }
  
  void Dynamics::updateAfter(Particle& particle)
  {
    //cout << "Dynamics::updateAfter" <<endl;
    particle.momentum() += timeStep_ * particle.force() / 2;
    particle.kineticEnergy() = pow(particle.momentum().norm(), 2) / particle.mass() / 2;
  }
  
  void Dynamics::updateAllControlVariates(Output& output, vector<Particle> const& configuration, size_t indexOfIteration) const
  {
    dvec q = configuration[0].position();
    dvec p = configuration[0].momentum();
    //assert(configuration.size() == 1);
    output.velocityCV()->update(p(0), generatorOn(output.velocityCV(), q, p), q, p, indexOfIteration);
    output.forceCV()->update(potential_->derivative(q)(0), generatorOn(output.forceCV(), q, p), q, p, indexOfIteration);
    output.lengthCV()->update(q(0), generatorOn(output.lengthCV(), q, p), q, p, indexOfIteration);
  }
  
  
  //#### Hamiltonian ####
  
  Hamiltonian::Hamiltonian(Input const& input, int const& indexOfReplica):Dynamics(input, indexOfReplica)
  {
  }
  
  double Hamiltonian::generatorOn(ControlVariate const* controlVariate, dvec const& position, dvec const& momentum) const
  {
    return momentum.dot(controlVariate->gradientQ(position, momentum))
    + force(position).dot(controlVariate->gradientP(position, momentum));   
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
  
  double Langevin::generatorOn(ControlVariate const* controlVariate, dvec const& position, dvec const& momentum) const
  {
    return momentum.dot(controlVariate->gradientQ(position, momentum))
    + force(position).dot(controlVariate->gradientP(position, momentum))
    + gamma_ * (- momentum.dot(controlVariate->gradientP(position, momentum)) 
		+ controlVariate->laplacienP(position, momentum) / beta_ );    
  }
  
    
  //#### Overdamped ####

  Overdamped::Overdamped(Input const& input, int const& indexOfReplica):
    UniformStochasticDynamics(input, indexOfReplica)
  {}
  
    
  void Overdamped::updateBefore(Particle& particle)
  {
    
  }
  
  // /!\ a changer dans le cas N != 1
  void Overdamped::updateAfter(Particle& particle)
  {
    //particle.position() += timeStep_ * particle.force() + sqrt(2*timeStep_/beta_) * rng_->gaussian();
  
    //assert(particle.force()(0) == force(particle.position())(0));
    dvec randomTerm = sqrt(2*timeStep_/beta_) * rng_->gaussian();
    dvec qtilde = particle.position() + .5 * timeStep_ * particle.force() + .5 * randomTerm;
    particle.position() += timeStep_ * force(qtilde) + randomTerm;
  }
  
  double Overdamped::generatorOn(ControlVariate const* controlVariate, dvec const& position, dvec const& momentum) const
  {
    return controlVariate->laplacienQ(position, momentum) / beta_
    + force(position).dot(controlVariate->gradientQ(position, momentum));    
  }

  
  //#### BoundaryStochasticDynamics ####

  BoundaryStochasticDynamics::BoundaryStochasticDynamics(Input const& input, int const& indexOfReplica):
    StochasticDynamics(input, indexOfReplica),
    betaLeft_(input.betaLeft()),
    betaRight_(input.betaRight()),
    temperatureLeft_(1/betaLeft_),
    temperatureRight_(1/betaRight())
  {} 
  
    double const& BoundaryStochasticDynamics::betaLeft() const
  {
    return betaLeft_;
  }
  
  double const& BoundaryStochasticDynamics::betaRight() const
  {
    return betaRight_;
  }
  
    double const& BoundaryStochasticDynamics::temperatureLeft() const
  {
    return temperatureLeft_;
  }
  
  double const& BoundaryStochasticDynamics::temperatureRight() const
  {
    return temperatureRight_;
  }

  //#### BoundaryLangevin ####
  
  BoundaryLangevin::BoundaryLangevin(const Input& input, const int& indexOfReplica):
    BoundaryStochasticDynamics(input, indexOfReplica),
    gamma_(input.gamma())
  {}
    
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
  
  //##BoundaryOverdamped ####
  
  BoundaryOverdamped::BoundaryOverdamped(const Input& input, const int& indexOfReplica):BoundaryStochasticDynamics(input, indexOfReplica){}
  
  /*void BoundaryOverdamped::updateBefore(Particle& particle)
  {
  }*/
  
  void BoundaryOverdamped::updateAfterLeft(Particle& particle)
  {
    particle.position() += timeStep_ * particle.force() + sqrt(2*timeStep_/betaLeft_) * rng_->gaussian();
  }
  
  void BoundaryOverdamped::updateAfterRight(Particle& particle)
  {
    particle.position() += timeStep_ * particle.force() + sqrt(2*timeStep_/betaRight_) * rng_->gaussian();
  }


    
}

#endif
