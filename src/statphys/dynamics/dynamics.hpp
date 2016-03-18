#ifndef SIMOL_DYNAMICS_HPP
#define SIMOL_DYNAMICS_HPP

#include "potential.hpp"
#include "particle.hpp"
#include "input.hpp"
#include "output.hpp"
#include "core/random/RNG.hpp"
#include "controlVariate.hpp"
# include <iostream>
#include "galerkin.hpp"

namespace simol
{
  class Dynamics;

  Dynamics* createDynamics(Input  const& input);

  class Dynamics
  {
    friend Dynamics* createDynamics(Input  const& input);
    public:
      Dynamics(Input const&  input);

			// Accessors
      double& timeStep();
      const double& timeStep() const;
      size_t& nbOfIterations();
      const size_t& nbOfIterations() const;
      double finalTime() const;
			size_t& nbOfThermalIterations();
      const size_t& nbOfThermalIterations() const;
			size_t& nbOfBurningIterations();
      const size_t& nbOfBurningIterations() const;
      const std::shared_ptr<RNG> rng() const;
      std::shared_ptr<RNG> rng();
      
      Vector<double>& externalForce() ;
      const Vector<double>& externalForce() const;
			double& externalForce(const int& i);
      const double& externalForce(const int& i) const;
			Galerkin* galerkin();
      
      //string name() {cout << "Dynamics" << endl;}

			// Accessors for daughters, fails if inadequate

			virtual const double& gamma() const {assert(false);}
			virtual double temperature() const {assert(false);}
      virtual const double& temperatureLeft() const {assert(false);}
      virtual const double& temperatureRight() const {assert(false);}
      virtual const double& betaLeft() const {assert(false);}
      virtual const double& betaRight() const {assert(false);}
      virtual double deltaTemperature() const {assert(false);}
      virtual const double& tauBending() const {assert(false);}
      virtual const double& xi() const {assert(false);}
			virtual double& xi() {assert(false);}
			int xiNbOfIterations() {assert(false);}
			
			virtual void initializeCountdown(Particle& /*particle*/){assert(false);};

      void resetForce(Particle& particle) const;
      
      virtual void updateBefore(Particle& particle);
      virtual void updateAfter(Particle& particle);
			virtual void updateOrsteinUhlenbeck(Particle& particle, double localBeta);
			virtual bool doMomentaExchange() const {return false;};
			virtual void updateMomentaExchange(Particle& /*particle1*/, Particle& /*particle2*/){assert(false);};
      virtual void bending(Particle& /*particle1*/, Particle& /*particle2*/) const {};
	protected:
      double timeStep_;
      size_t nbOfIterations_, nbOfThermalIterations_, nbOfBurningIterations_;
      //double timeStep;
      Vector<double> externalForce_;
			std::shared_ptr<RNG> rng_;
			Galerkin* galerkin_;
  };

  class Hamiltonian : public Dynamics
  {
  public:
    Hamiltonian(Input const&  input);
  };

  class StochasticDynamics : public Dynamics
  {
		double xi_;
  public:
    StochasticDynamics(Input const& input);
		virtual const double& xi() const;
		virtual double& xi();
		int xiNbOfIterations();
		virtual bool doMomentaExchange() const;
		virtual void initializeCountdown(Particle& particle);
		virtual void updateMomentaExchange(Particle& particle1, Particle& particle2);
  };

  class UniformStochasticDynamics : public StochasticDynamics
  {
  protected:
    double beta_;
    double temperature_;
  public:
    UniformStochasticDynamics(Input const& input);
    void initializeMomenta(vector<Particle>& configuration);
    virtual double  temperature() const;
		virtual const double& temperatureLeft() const;
    virtual const double& temperatureRight() const;
    virtual const double& beta() const;
		virtual const double& betaLeft() const;
    virtual const double& betaRight() const;
  };

  class Langevin : public UniformStochasticDynamics
  {
  public:
    Langevin(Input const& input);
    virtual const double& gamma() const;
    double sigma() const;
    virtual void updateAfter(Particle& particle);
  protected:
    double gamma_;
  };

  class Overdamped : public UniformStochasticDynamics
  {
  public:
    Overdamped(Input const& input);
    virtual void updateBefore(Particle& particle);
    virtual void updateAfter(Particle& particle);
 };

  class BoundaryLangevin : public StochasticDynamics
  {
  public:
    BoundaryLangevin(Input const& input);
    virtual const double& betaLeft() const;
    virtual const double& betaRight() const;
		virtual double temperature() const;
    virtual const double& temperatureLeft() const;
    virtual const double& temperatureRight() const;
    double deltaTemperature() const;
		const double& gamma() const;
    double sigmaLeft() const;
    double sigmaRight() const;
    const double& tauBending() const;

    void initializeMomenta(vector<Particle>& configuration);
		virtual void bending(Particle& particle1, Particle& particle2) const;
  protected:
    double betaLeft_;
    double betaRight_;
    double temperatureLeft_;
    double temperatureRight_;
		double gamma_;
    double tauBending_;
  };


}

#endif
