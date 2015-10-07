#ifndef SIMOL_DYNAMICS_HPP
#define SIMOL_DYNAMICS_HPP

#include "potential.hpp"
#include "particle.hpp"
#include "input.hpp"
#include "RNG.hpp"
# include <iostream>

namespace simol
{
  class Dynamics;
  
  Dynamics* createDynamics(Input  const& input);
    
  class Dynamics
  {
    public:
      Dynamics(Input const&  input);
      virtual ~Dynamics();

      Potential const & potential() const;
      double potential(dvec const& position) const;
      double potential(double const& position) const;
      dvec force(dvec const& position) const;
      //friend Dynamics* createDynamics(Potential const& potential);
      friend Dynamics* createDynamics(Input  const& input);
      virtual void update(Particle& particle,  double const timeStep) = 0;
      void resetForce(Particle& particle) const;
      void computeForce(Particle& particle) const;
      void interaction(Particle& particle1, Particle& particle2) const;
    private:
      Potential* potential_;
      //double timeStep;
  };
  
  class Hamiltonian : public Dynamics
  {
  public:
    Hamiltonian(Input const&  input);
    void update(Particle& particle,  double const timeStep);
  };
  
  
  class Langevin : public Dynamics
  {
  public:
    Langevin(Input const& input);
    double const& temperature() const;
    double const& beta() const;
    double const& gamma() const;
    double const& sigma() const;
    void update(Particle& particle,  double const timeStep);
  private:
    double temperature_;
    double beta_;
    double gamma_;
    double sigma_;
    RNG rng_;
  };
  
  class Overdamped : public Dynamics
  {
  public:
    Overdamped(Input const& input);
    double const& temperature() const;
    double const& beta() const;
    void update(Particle& particle,  double const timeStep);
  private:
    double temperature_;
    double beta_;
    RNG rng_;
  };
  


}

//#include "dynamics.cpp"

#endif
