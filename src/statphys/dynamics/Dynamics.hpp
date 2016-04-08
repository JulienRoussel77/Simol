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
    
    // Accessors
    double& timeStep();
    const double& timeStep() const;
    size_t& nbOfIterations();
    size_t const& nbOfIterations() const;
    double finalTime() const;
    size_t& nbOfThermalIterations();
    size_t const& nbOfThermalIterations() const;
    size_t& nbOfBurningIterations();
    size_t const& nbOfBurningIterations() const;
    std::shared_ptr<RNG> const& rng() const;
    std::shared_ptr<RNG>& rng();
    virtual const double&  temperature() const;
    virtual const double& temperatureLeft() const;
    virtual const double& temperatureRight() const;
    virtual double deltaTemperature() const;
    virtual const double& beta() const;
    virtual const double& betaLeft() const;
    virtual const double& betaRight() const;
    //Vector<double>& externalForce() ;
    //Vector<double> const& externalForce() const;
    double& externalForce(const int& i);
    double const& externalForce(const int& i) const;
    Galerkin* galerkin();
    
      
    //void resetForce(Particle& particle) const;
    
    void verletFirstPart(Particle& particle);
    void verletSecondPart(Particle& particle);
    virtual void updateBefore(Particle& particle);
    virtual void updateAfter(Particle& particle);
    virtual bool doMomentaExchange() const {return false;};
    virtual void updateMomentaExchange(Particle& /*particle1*/, Particle& /*particle2*/){assert(false);};
    virtual void bending(Particle& /*particle1*/, Particle& /*particle2*/) const {};
  protected:
    double timeStep_;
    size_t nbOfIterations_, nbOfThermalIterations_, nbOfBurningIterations_;
    double beta_;
    double temperature_;
    //double timeStep;
    std::shared_ptr<RNG> rng_;
    Galerkin* galerkin_;
  };

}

#endif
