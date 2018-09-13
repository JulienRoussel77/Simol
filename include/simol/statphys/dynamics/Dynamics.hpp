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
#include "simol/statphys/controlVariate/Operator.hpp"

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
      double const& timeStep() const;
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
      virtual double const&  temperature() const;
      virtual double const& temperatureLeft() const;
      virtual double const& temperatureRight() const;
      virtual double deltaTemperature() const;
      virtual double const& beta() const;
      virtual double const& betaLeft() const;
      virtual double const& betaRight() const;
      //virtual bool const& isMollified() const;
      
      virtual void sampleInternalEnergies(System&) const {};
      virtual void computeOutput(System const& syst, Output& output, long int iOfStep) const;
      virtual void computeControlVariate(System const& syst, Output& output) const;
      virtual void writeOutput(System const& syst, Output& output, long int iOfStep) const;
      virtual void writeFinalOutput(System const& syst, Output& output) const;
      virtual void thermalize(System& syst) const;
      virtual void sampleSystem(System& syst) const;
      virtual void launch(System& syst, Output& output);
      
      virtual void simulate(System& /*syst*/) const {cout << "simulate(System&)" << endl;};
       
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
      
      
      
      /*virtual double& lagrangeMultiplier() {throw runtime_error("resetConstraint not implemented for this system !");}
      virtual double const& lagrangeMultiplier() const {throw runtime_error("resetConstraint not implemented for this system !");}*/
      
      //-- output functions --
      virtual void computeKineticEnergy(Output& output, System const& syst) const;
      virtual void computePotentialEnergy(Output& output, System const& syst) const;
      virtual void computePressure(Output& output, System const& syst) const;
      virtual void computeInternalEnergy(Output& output, System const& syst) const;
      virtual void computeInternalTemperature(Output& output, System const& syst) const;
      
      virtual void computeProfileBiChain(Output&, System const&, long int) const {}
      virtual void computeProfileTriChain(Output&, System const&, long int) const {}
            
      //--for output --
      virtual void getThermo(Output& output) const;
      virtual void getPressure(Output& output) const;
    protected:
      // contains all the parameters of the dynamics
      DynamicsParameters parameters_;
      // timestep of the integrator
      double timeStep_;
      // number of steps of integration with output
      long int nbOfSteps_;
      // number of steps of thermalization (see method thermalize)
      long int thermalizationNbOfSteps_;
      // number of steps of integration without output
      long int burninNbOfSteps_;
      // random number generator
      std::shared_ptr<RNG> rng_;
  };

}

#endif
