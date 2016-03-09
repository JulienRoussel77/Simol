#ifndef SIMOL_DYNAMICS_HPP
#define SIMOL_DYNAMICS_HPP

#include "potential.hpp"
#include "particle.hpp"
#include "input.hpp"
#include "output.hpp"
#include "core/random/RNG.hpp"
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
      double const& timeStep() const;
      size_t& numberOfIterations();
      const size_t& numberOfIterations() const;
      double finalTime() const;
      Potential& potential();
      double potential(Vector<double> const& position) const;
      double potential(double const& position) const;
      Vector<double> force(Vector<double> const& position) const;
      Vector<double>& externalForce() ;
      const Vector<double>& externalForce() const;
      double& externalForce(int const& i);
      double const& externalForce(int const& i) const;
			virtual double const& gamma() const {assert(false);}
			virtual double temperature() const {assert(false);}
      virtual double const& temperatureLeft() const {assert(false);}
      virtual double const& temperatureRight() const {assert(false);}
      virtual double deltaTemperature() const {assert(false);}
      virtual double const& tauBending() const {assert(false);}
      //friend Dynamics* createDynamics(Potential const& potential);

      virtual void initializeMomenta(vector<Particle>& configuration);
			virtual Vector<double> drawMomentum(double /*localBeta*/, double /*mass*/){assert(false);return Vector<double>();};
      virtual void setRNG(RNG* /*rng*/){};
      void resetForce(Particle& particle) const;
      void computeForce(Particle& particle) const;
      void interaction(Particle& particle1, Particle& particle2) const;
      void triInteraction(Particle& particle1, Particle& particle2, Particle& particle3) const;
      virtual void updateBefore(Particle& particle);
      virtual void updateAfter(Particle& particle);
      virtual void updateAfterLeft(Particle& /*particle*/) {};
      virtual void updateAfterRight(Particle& /*particle*/) {};
      virtual void bending(Particle& /*particle1*/, Particle& /*particle2*/) const {};
      virtual MatrixXd generatorOn(ControlVariate const* /*controlVariate*/, vector<Particle> const& /*configuration*/) const{cout << "operator not implemented !"; assert(false); return MatrixXd();}
      virtual void updateAllControlVariates(Output& output, vector<Particle> const& configuration, size_t indexOfIteration) const;
			virtual double computeMeanPotLaw(double betaLocal) const;
	protected:
      double timeStep_;
      size_t numberOfIterations_;
      Potential* potential_;
      //double timeStep;
      Vector<double> externalForce_;
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
		virtual Vector<double> drawMomentum(double localBeta, double mass);
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
    virtual double  temperature() const;
    virtual double const& beta() const;
    virtual double const& temperatureLeft() const;
    virtual double const& temperatureRight() const;
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
  };

  class BoundaryStochasticDynamics : public StochasticDynamics
  {
  public:
    BoundaryStochasticDynamics(Input const& input, int const& indexOfReplica=1);
    virtual double const& betaLeft() const;
    virtual double const& betaRight() const;
		virtual double temperature() const;
    virtual double const& temperatureLeft() const;
    virtual double const& temperatureRight() const;
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

    double const& tauBending() const;

    //virtual void updateAfter(Particle& particle);
    virtual void updateAfterLeft(Particle& particle);
    virtual void updateAfterRight(Particle& particle);
    virtual void bending(Particle& particle1, Particle& particle2) const;
    MatrixXd generatorOn(ControlVariate const* controlVariate, vector<Particle> const& configuration) const;
  protected:
    double gamma_;
    double tauBending_;
  };


}

//#include "dynamics.cpp"

#endif
