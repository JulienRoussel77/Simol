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
    particle.kineticEnergy() = pow(particle.momentum().norm(), 2) / particle.mass() / 2;
    particle.potentialEnergy() = potential(particle.position());
    particle.force() = force(particle.position());
  }
  
  void Dynamics::interaction(Particle& particle1, Particle& particle2) const
  {
    dvec r12 = particle2.position() - particle1.position();
    double energy12 = potential(r12);
    dvec force12 = force(r12);
    
    particle1.potentialEnergy() += energy12;
    particle1.potentialEnergy() += energy12;
    particle1.force() += force12;
    particle1.force() -= force12;
  }
  
  Hamiltonian::Hamiltonian(Input const& input, int const& indexOfReplica):Dynamics(input, indexOfReplica)
  {
  }

  void Hamiltonian::update(Particle& particle)
  {
    verlet_scheme(particle, potential(), timeStep());
  }

  //#### Langevin ####

  Langevin::Langevin(Input const& input, int const& indexOfReplica):
    Dynamics(input, indexOfReplica), 
    temperature_(input.temperature()), 
    beta_(1/temperature_), 
    gamma_(input.gamma()), 
    sigma_(2*gamma_/beta_)
  {}
  

  
    double const& Langevin::temperature() const
    {
      return temperature_;
    }
      
    double const& Langevin::beta() const
    {
      return beta_;
    }
    
    double const& Langevin::gamma() const
    {
      return gamma_;
    }
    
    double const& Langevin::sigma() const
    {
      return sigma_;
    }
    
    void Langevin::setRNG(RNG* rng)
    {
      rng_ = rng;
    };
    
    void Langevin::update(Particle& particle)
    {
      verlet_scheme(particle, potential(), timeStep());
      exact_OU_scheme(particle, gamma_, beta_, timeStep(), rng_->gaussian());
    }
    
  //#### Overdamped ####

  Overdamped::Overdamped(Input const& input, int const& indexOfReplica):
    Dynamics(input, indexOfReplica), 
    temperature_(input.temperature()), 
    beta_(1/temperature_)
  {}
  

  
    double const& Overdamped::temperature() const
    {
      return temperature_;
    }
      
    double const& Overdamped::beta() const
    {
      return beta_;
    }
    
    void Overdamped::setRNG(RNG* rng)
    {
      rng_ = rng;
    };
    
    void Overdamped::update(Particle& particle)
    {
      maruyama_scheme(particle, beta_, potential(), timeStep(), rng_->gaussian());
    }
    
}

#endif
