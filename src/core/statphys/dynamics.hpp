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

      Potential const & potential() const;
      //friend Dynamics* createDynamics(Potential const& potential);
      friend Dynamics* createDynamics(Input  const& input);
      virtual void update(Particle& particle,  double const timeStep) = 0;
    private:
      Potential potential_;

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
  


}

//#include "dynamics.cpp"

#endif
