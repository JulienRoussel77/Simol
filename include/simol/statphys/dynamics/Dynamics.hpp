#ifndef SIMOL_DYNAMICS_HPP
#define SIMOL_DYNAMICS_HPP

#include "simol/statphys/potential/Potential.hpp"
#include "simol/statphys/system/Particle.hpp"
#include "simol/statphys/input/Input.hpp"
#include "simol/statphys/output/Output.hpp"
#include "simol/core/random/RNG.hpp"
#include "simol/statphys/controlVariate/ControlVariate.hpp"
#include "simol/statphys/controlVariate/CVBasis.hpp"
# include <iostream>
#include "simol/statphys/controlVariate/AllGalerkins.hpp"

#include "simol/statphys/system/System.hpp"
#include "simol/statphys/dynamics/DynamicsParameters.hpp"

namespace simol
{

  class Dynamics
  {

    public:
      Dynamics(Input const&  input);
      virtual string dynamicsName() const=0;


      //Contains all the parameters of the dynamics
      DynamicsParameters const& parameters() const {return parameters_;}
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
      
      //--for output --
      virtual void getThermo(Output& output) const = 0;
      virtual void getPressure(Output& output) const = 0;
       
      //-- for DPDE --
      virtual double internalTemperature(double /*intEnergy*/) const {return 0;}
      virtual double entropy_derivative(double /*intEnergy*/) const {return 0;}
      virtual void specificComputeOutput(Output& /*output*/) const {}

    //-- external forces --
      double& externalForce(const int& i);
      double const& externalForce(const int& i) const;

      //-- for numerical integration --
      void verletFirstPart(Particle& particle) const;
      void verletSecondPart(Particle& particle) const;
      void updateMomentum(Particle& particle) const;
      void updatePosition(Particle& particle) const;
      virtual bool doMomentaExchange() const {return false;};
      virtual void updateMomentaExchange(Particle& /*particle1*/, Particle& /*particle2*/) {assert(false);};
      virtual void bending(Particle& /*particle1*/, Particle& /*particle2*/) const {};
      
      //-- control variates --
      //Galerkin* galerkin();
      //shared_ptr<CVBasis> createCvBasis(Input const& input);
      //shared_ptr<CVBasis> cvBasis();
      //return CVBasis(dynamic_cast<TensorBasis*>(&basis_), make_shared<DVec>(CVcoeffsVec()));
      
      virtual void computeGeneratorOnBasis(shared_ptr<CVBasis>, System const&) const {};
      
      
      virtual double& lagrangeMultiplier() {throw runtime_error("resetConstraint not implemented for this system !");}
      virtual const double& lagrangeMultiplier() const {throw runtime_error("resetConstraint not implemented for this system !");}
      
      //-- output functions --
      virtual void computeKineticEnergy(Output& output, System const& syst) const;
      virtual void computePotentialEnergy(Output& output, System const& syst) const;
      virtual void computePressure(Output& output, System const& syst) const;
      virtual void computeInternalEnergy(Output& output, System const& syst) const;
      virtual void computeInternalTemperature(Output& output, System const& syst) const;
      
      virtual void computeProfileBiChain(Output&, System const&, long int) const {}
      virtual void computeProfileTriChain(Output&, System const&, long int) const {}
    protected:
      DynamicsParameters parameters_;
      double timeStep_;
      long int nbOfSteps_, thermalizationNbOfSteps_, burninNbOfSteps_;
      double beta_;
      double temperature_;
      std::shared_ptr<RNG> rng_;
  };

}

#endif
