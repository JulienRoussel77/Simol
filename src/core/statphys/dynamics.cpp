#ifndef SIMOL_DYNAMICS_IPP
#define SIMOL_DYNAMICS_IPP

#include "dynamics.hpp"


namespace simol
{

    Dynamics* createDynamics(Input const& input)
  {
    if (input.dynamicsName() == "Hamiltonian")
      return new Hamiltonian(input);
    else if (input.dynamicsName() == "Langevin")
      return new Langevin(input);
    else if (input.dynamicsName() == "Overdamped")
      return new Overdamped(input);
    else
      std::cout << input.dynamicsName() << " is not a valid dynamics !" << std::endl;
    
    return 0;
  }
  
  Dynamics::Dynamics(Input const& input)
  {
    potential_ = createPotential(input);
  }
  
  Dynamics::~Dynamics()
  {
    delete potential_;
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
   return potential_->force(position); 
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
  
  Hamiltonian::Hamiltonian(Input const& input):Dynamics(input)
  {
  }

  void Hamiltonian::update(Particle& particle,  double const timeStep)
  {
    verlet_scheme(particle, potential(), timeStep);
  }

  //#### Langevin ####

  Langevin::Langevin(Input const& input):
    Dynamics(input), 
    temperature_(input.temperature()), 
    beta_(1/temperature_), 
    gamma_(input.gamma()), 
    sigma_(2*gamma_/beta_), 
    rng_(1, input.dimension())
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
    
    void Langevin::update(Particle& particle,  double const timeStep)
    {
      verlet_scheme(particle, potential(), timeStep);
      exact_OU_scheme(particle, gamma_, beta_, timeStep, rng_.gaussian());
    }
    
  //#### Overdamped ####

  Overdamped::Overdamped(Input const& input):
    Dynamics(input), 
    temperature_(input.temperature()), 
    beta_(1/temperature_), 
    rng_(1, input.dimension())
  {}
  

  
    double const& Overdamped::temperature() const
    {
      return temperature_;
    }
      
    double const& Overdamped::beta() const
    {
      return beta_;
    }
    
    void Overdamped::update(Particle& particle,  double const timeStep)
    {
      maruyama_scheme(particle, beta_, potential(), timeStep, rng_.gaussian());
    }
    
}

#endif
