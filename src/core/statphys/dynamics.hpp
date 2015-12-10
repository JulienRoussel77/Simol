#ifndef SIMOL_DYNAMICS_HPP
#define SIMOL_DYNAMICS_HPP

#include "potential.hpp"
#include "particle.hpp"
#include "input.hpp"
#include "output.hpp"
#include "RNG.hpp"
#include "controlVariate.hpp"
# include <iostream>

namespace simol
{
  class Dynamics;
  
  Dynamics* createDynamics(Input  const& input, size_t indexOfReplica=1);
    
  class Dynamics
  {
    friend Dynamics* createDynamics(Input  const& input, size_t indexOfReplica);
    public:
      Dynamics(Input const&  input, int const& indexOfReplica=1);
      virtual ~Dynamics();

      Potential const & potential() const;
      double& timeStep();
      const double& timeStep() const;
      size_t& numberOfIterations();
      const size_t& numberOfIterations() const;
      double finalTime() const;
      Potential& potential();
      double potential(dvec const& position) const;
      double potential(double const& position) const;
      dvec force(dvec const& position) const;
      dvec& externalForce() ;
      const dvec& externalForce() const;
      double& externalForce(int const& i) ;
      const double& externalForce(int const& i) const;
      virtual const double& temperatureLeft() const {assert(false);}
      virtual const double& temperatureRight() const {assert(false);}
      //friend Dynamics* createDynamics(Potential const& potential);
      
      virtual void initializeMomenta(vector<Particle>& configuration);
      virtual void setRNG(RNG* rng){};
      void resetForce(Particle& particle) const;
      void computeForce(Particle& particle) const;
      void interaction(Particle& particle1, Particle& particle2) const;
      void triInteraction(Particle& particle1, Particle& particle2, Particle& particle3) const;
      virtual void updateBefore(Particle& particle);
      virtual void updateAfter(Particle& particle);
      virtual void updateAfterLeft(Particle& particle) {};
      virtual void updateAfterRight(Particle& particle) {};
      virtual MatrixXd generatorOn(ControlVariate const* controlVariate, vector<Particle> const& configuration) const{cout << "operator not implemented !"; assert(false); return MatrixXd();}
      virtual void updateAllControlVariates(Output& output, vector<Particle> const& configuration, size_t indexOfIteration) const;
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
    virtual MatrixXd generatorOn(ControlVariate const* controlVariate, vector<Particle> const& configuration) const;

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
  protected:
    double beta_;
    double temperature_;
  public:
    UniformStochasticDynamics(Input const& input, int const& indexOfReplica=1);
    void initializeMomenta(vector<Particle>& configuration);
    virtual double const& temperature() const;
    virtual double const& beta() const;
    virtual const double& temperatureLeft() const;
    virtual const double& temperatureRight() const;
  };
  
  class Langevin : public UniformStochasticDynamics
  {
  public:
    Langevin(Input const& input, int const& indexOfReplica=1);
    double const& gamma() const;
    double sigma() const;
    virtual void updateAfter(Particle& particle);
    virtual MatrixXd generatorOn(ControlVariate const* controlVariate, vector<Particle> const& configuration) const;
  protected:
    double gamma_;
  };
  
  class Overdamped : public UniformStochasticDynamics
  {
  public:
    Overdamped(Input const& input, int const& indexOfReplica=1);
    virtual void updateBefore(Particle& particle);
    virtual void updateAfter(Particle& particle);
    virtual MatrixXd generatorOn(ControlVariate const* controlVariate, vector<Particle> const& configuration) const;
    //virtual std::function<double (ControlVariate const*, dvec const&, dvec const&)> generator() const;
  };
  
  class BoundaryStochasticDynamics : public StochasticDynamics
  {
  public:
    BoundaryStochasticDynamics(Input const& input, int const& indexOfReplica=1);
    virtual const double& betaLeft() const;
    virtual const double& betaRight() const;
    virtual const double& temperatureLeft() const;
    virtual const double& temperatureRight() const;
    double deltaTemperature() const;
    
    //virtual void updateAfter(Particle& particle);
    virtual void updateAfterLeft(Particle& particle) = 0;
    virtual void updateAfterRight(Particle& particle) = 0;
    
    void initializeMomenta(vector<Particle>& configuration);
  protected:
    double betaLeft_; 
    double betaRight_; 
    double temperatureLeft_;
    double temperatureRight_;
  };  
  
  
  class BoundaryLangevin : public BoundaryStochasticDynamics
  {
  public:
    BoundaryLangevin(Input const& input, int const& indexOfReplica=1);
    double const& gamma() const;
    double sigmaLeft() const;
    double sigmaRight() const;
    
    //virtual void updateAfter(Particle& particle);
    virtual void updateAfterLeft(Particle& particle);
    virtual void updateAfterRight(Particle& particle);
    MatrixXd generatorOn(ControlVariate const* controlVariate, vector<Particle> const& configuration) const;
  protected:
    double gamma_;
  };


}

//#include "dynamics.cpp"

#endif