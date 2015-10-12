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
  
  Dynamics* createDynamics(Input  const& input, int const& indexOfReplica=1);
    
  class Dynamics
  {
    friend Dynamics* createDynamics(Input  const& input, int const& indexOfReplica);
    public:
      Dynamics(Input const&  input, int const& indexOfReplica=1);
      virtual ~Dynamics();

      Potential const & potential() const;
      double& timeStep();
      const double& timeStep() const;
      size_t& numberOfIterations();
      const size_t& numberOfIterations() const;
      double finalTime() const;
      double potential(dvec const& position) const;
      double potential(double const& position) const;
      dvec force(dvec const& position) const;
      dvec& externalForce() ;
      const dvec& externalForce() const;
      double& externalForce(int const& i) ;
      const double& externalForce(int const& i) const;
      //friend Dynamics* createDynamics(Potential const& potential);
      virtual void setRNG(RNG* rng){};
      virtual void update(Particle& particle) = 0;
      void resetForce(Particle& particle) const;
      void computeForce(Particle& particle) const;
      void interaction(Particle& particle1, Particle& particle2) const;
    private:
      double timeStep_;
      size_t numberOfIterations_;
      Potential* potential_;
      //double timeStep;
      dvec externalForce_;
  };
  
  class Hamiltonian : public Dynamics
  {
  public:
    Hamiltonian(Input const&  input, int const& indexOfReplica=1);
    void update(Particle& particle);
  };
  
  
  class Langevin : public Dynamics
  {
  public:
    Langevin(Input const& input, int const& indexOfReplica=1);
    double const& temperature() const;
    double const& beta() const;
    double const& gamma() const;
    double const& sigma() const;
    void setRNG(RNG* rng);
    void update(Particle& particle);
  private:
    double temperature_;
    double beta_;
    double gamma_;
    double sigma_;
    RNG* rng_;
  };
  
  class Overdamped : public Dynamics
  {
  public:
    Overdamped(Input const& input, int const& indexOfReplica=1);
    double const& temperature() const;
    double const& beta() const;
    void setRNG(RNG* rng);
    void update(Particle& particle);
  private:
    double temperature_;
    double beta_;
    RNG* rng_;
  };
  


}

//#include "dynamics.cpp"

#endif
