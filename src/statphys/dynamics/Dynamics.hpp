#ifndef SIMOL_DYNAMICS_HPP
#define SIMOL_DYNAMICS_HPP

#include "Potential.hpp"
#include "Particle.hpp"
#include "Input.hpp"
#include "Output.hpp"
#include "core/random/RNG.hpp"
#include "ControlVariate.hpp"
# include <iostream>
#include "Galerkin.hpp"

namespace simol
{

  class Dynamics
  {
    
  public:
    Dynamics(Input const&  input);
    virtual ~Dynamics() = default;
    
    virtual void printName() const;
    
    //-- time steps and iterations --
    double& timeStep();
    const double& timeStep() const;
    int& nbOfIterations();
    int const& nbOfIterations() const;
    double finalTime() const;
    int& nbOfThermalIterations();
    int const& nbOfThermalIterations() const;
    int& nbOfBurnInIterations();
    int const& nbOfBurnInIterations() const;
    
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
    virtual void updateMomentaExchange(Particle& /*particle1*/, Particle& /*particle2*/){assert(false);};
    virtual void bending(Particle& /*particle1*/, Particle& /*particle2*/) const {};
    
    //-- control variates --
    Galerkin* galerkin();
    
  protected:
    
    double timeStep_;
    int nbOfIterations_, nbOfThermalIterations_, nbOfBurnInIterations_;
    double beta_;
    double temperature_;
    std::shared_ptr<RNG> rng_;
    Galerkin* galerkin_;
  };

}

#endif
