#ifndef SIMOL_DYNAMICS_IPP
#define SIMOL_DYNAMICS_IPP

#include "dynamics.hpp"


namespace simol
{

  Dynamics::Dynamics(Input const& input) : 
    potential_(input.potParameter(), 2*M_PI/input.length())
  {}
  
  /*Dynamics* createDynamics(Potential const& potential)
  {
    return new Hamiltonian(potential);
  }*/
  
  Dynamics* createDynamics(Input const& input)
  {
    if (input.name() == "Hamiltonian")
      return new Hamiltonian(input);
    else if (input.name() == "Langevin")
      return new Langevin(input);
    else
      std::cout << "Method not valid !" << std::endl;
    
    return 0;
  }

  Potential const & Dynamics::potential() const
  { return potential_; }
  
  void Dynamics::resetForce(Particle& particle) const
  {
    particle.potentialEnergy() = 0;
    particle.force() = dvec(3, 0);
  }

  void Dynamics::computeForce(Particle& particle) const
  {
    particle.kineticEnergy() = pow(particle.momentum().norm(), 2) / particle.mass() / 2;
    particle.potentialEnergy() = potential_(particle.position());
    particle.force() = potential_.force(particle.position());
  }
  
  Hamiltonian::Hamiltonian(Input const& input):Dynamics(input)
  {}

  void Hamiltonian::update(Particle& particle,  double const timeStep)
  {
    verlet_scheme(particle, potential(), timeStep);
  }


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
}

#endif
