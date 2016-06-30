#ifndef SIMOL_DYNAMICS_HPP
#define SIMOL_DYNAMICS_HPP

#include "simol/statphys/potential/Potential.hpp"
#include "simol/statphys/system/Particle.hpp"
#include "simol/statphys/input/Input.hpp"
#include "simol/statphys/output/Output.hpp"
#include "simol/core/random/RNG.hpp"
#include "simol/statphys/controlVariate/ControlVariate.hpp"
# include <iostream>
#include "simol/statphys/controlVariate/Galerkin.hpp"

namespace simol
{

  class Dynamics
  {

    public:
      Dynamics(Input const&  input);
      virtual ~Dynamics() = default;

      virtual void printName() const;

      //-- time step and numbers of steps --
      double& timeStep();
      const double& timeStep() const;
      long int& nbOfSteps();
      long int const& nbOfSteps() const;
      double finalTime() const;
      long int& thermalizationNbOfSteps();
      long int const& thermalizationNbOfSteps() const;
      long int& burninNbOfSteps();
      long int const& burninNbOfSteps() const;

      //-- random numbers ---
      std::shared_ptr<RNG> const& rng() const;
      std::shared_ptr<RNG>& rng();

      //-- (inverse) temperatures --
      virtual const double&  temperature() const;
      virtual const double& temperatureLeft() const;
      virtual const double& temperatureRight() const;
      virtual double deltaTemperature() const;
      virtual const double& beta() const;
      virtual const double& betaLeft() const;
      virtual const double& betaRight() const;

      //-- external forces --
      double& externalForce(const int& i);
      double const& externalForce(const int& i) const;

      //-- for numerical integration --
      void verletFirstPart(Particle& particle);
      void verletSecondPart(Particle& particle);
      virtual void updateBefore(Particle& particle);
      virtual void updateAfter(Particle& particle);
      virtual bool doMomentaExchange() const {return false;};
      virtual void updateMomentaExchange(Particle& /*particle1*/, Particle& /*particle2*/) {assert(false);};
      virtual void bending(Particle& /*particle1*/, Particle& /*particle2*/) const {};

      //-- control variates --
      Galerkin* galerkin();

    protected:

      double timeStep_;
      long int nbOfSteps_, thermalizationNbOfSteps_, burninNbOfSteps_;
      double beta_;
      double temperature_;
      std::shared_ptr<RNG> rng_;
      Galerkin* galerkin_;
  };

}

#endif
