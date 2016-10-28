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
#include "simol/statphys/controlVariate/Galerkin.hpp"

#include "simol/statphys/system/System.hpp"
#include "simol/statphys/dynamics/DynamicsParameters.hpp"

namespace simol
{

  class Dynamics
  {

    public:
      Dynamics(Input const&  input);
      virtual ~Dynamics() {if (galerkin_) delete galerkin_;};
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
      void verletFirstPart(Particle& particle);
      void verletSecondPart(Particle& particle);
      //virtual void updateBefore(Particle& particle);
      //virtual void updateAfter(Particle& particle);
      virtual bool doMomentaExchange() const {return false;};
      virtual void updateMomentaExchange(Particle& /*particle1*/, Particle& /*particle2*/) {assert(false);};
      virtual void bending(Particle& /*particle1*/, Particle& /*particle2*/) const {};
      
      //-- control variates --
      Galerkin* galerkin();
      //shared_ptr<CVBasis> cvBasis() {return galerkin_?make_shared<CVBasis>(galerkin_->makeCvBasis()):nullptr;}
      shared_ptr<CVBasis> cvBasis() {return galerkin_?make_shared<CVBasis>(dynamic_cast<TensorBasis*>(&galerkin_->basis()), make_shared<DVec>(galerkin_->CVcoeffsVec())):nullptr;}
      //return CVBasis(dynamic_cast<TensorBasis*>(&basis_), make_shared<DVec>(CVcoeffsVec()));
      
      virtual void computeGeneratorOnBasis(CVBasis&, System const&) const {};
      
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
      Galerkin* galerkin_;
  };

}

#endif
