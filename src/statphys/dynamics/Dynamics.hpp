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

  class Dynamics
  {

    public:
      Dynamics(Input const&  input);
      virtual ~Dynamics() = default;

      virtual void printName() const;

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
      const std::shared_ptr<RNG>& rng() const;
      std::shared_ptr<RNG>& rng();

      Vector<double>& externalForce() ;
      const Vector<double>& externalForce() const;
	  double& externalForce(const int& i);
      const double& externalForce(const int& i) const;
	  Galerkin* galerkin();

	  // Accessors for daughters, fails if inadequate

	  virtual const double& gamma() const {assert(false);}
      virtual const double& temperatureLeft() const {assert(false);}
      virtual double deltaTemperature() const {assert(false);}

      void resetForce(Particle& particle) const;

      void verletFirstPart(Particle& particle);
      void verletSecondPart(Particle& particle);
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

}

#endif
