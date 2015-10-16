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
      void resetForce(Particle& particle) const;
      void computeForce(Particle& particle) const;
      void interaction(Particle& particle1, Particle& particle2) const;
      virtual void updateBefore(Particle& particle);
      virtual void updateAfter(Particle& particle);
      virtual void updateAfterLeft(Particle& particle) {updateAfter(particle);};
      virtual void updateAfterRight(Particle& particle) {updateAfter(particle);};
    protected:
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
  };
  
  class StochasticDynamics : public Dynamics
  {
  public:
    StochasticDynamics(Input const& input, int const& indexOfReplica=1);
    virtual void setRNG(RNG* rng);
  protected:
    RNG* rng_;
  };
  
  class UniformStochasticDynamics : public StochasticDynamics
  {
  public:
    UniformStochasticDynamics(Input const& input, int const& indexOfReplica=1);
    virtual double const& temperature() const;
    virtual double const& beta() const;
  protected:
    double temperature_;
    double beta_;
  };
  
  class Langevin : public UniformStochasticDynamics
  {
  public:
    Langevin(Input const& input, int const& indexOfReplica=1);
    double const& gamma() const;
    double sigma() const;
    virtual void updateAfter(Particle& particle);
  protected:
    double gamma_;
  };
  
  class Overdamped : public UniformStochasticDynamics
  {
  public:
    Overdamped(Input const& input, int const& indexOfReplica=1);
    virtual void updateBefore(Particle& particle);
    virtual void updateAfter(Particle& particle);
  };
  
  class BoundaryLangevin : public StochasticDynamics
  {
  public:
    BoundaryLangevin(Input const& input, int const& indexOfReplica=1);
    const double& betaLeft() const;
    const double& betaRight() const;
    const double& temperatureLeft() const;
    const double& temperatureRight() const;
    double const& gamma() const;
    double sigmaLeft() const;
    double sigmaRight() const;
    //virtual void updateAfter(Particle& particle);
    virtual void updateAfterLeft(Particle& particle);
    virtual void updateAfterRight(Particle& particle);
  protected:
    double betaLeft_; 
    double betaRight_; 
    double temperatureLeft_;
    double temperatureRight_;
    double gamma_;
  };
  


}

//#include "dynamics.cpp"

#endif
