#ifndef SIMOL_DYNAMICS_IPP
#define SIMOL_DYNAMICS_IPP

#include "simol/statphys/dynamics/Dynamics.hpp"

using std::cout;
using std::endl;

namespace simol
{



  ///
  ///Constructor
  Dynamics::Dynamics(Input const& input):
    timeStep_(input.timeStep()),
    nbOfSteps_(input.nbOfSteps()),
    thermalizationNbOfSteps_(input.thermalizationNbOfSteps()),
    burninNbOfSteps_(input.burninNbOfSteps()),
    beta_(input.beta()),
    temperature_(1/beta_),
    galerkin_(nullptr)
  {}

  void Dynamics::printName() const
  {
    std::cout << "DynamicsType = Dynamics" << std::endl;
  }




  ///
  ///Read-write accessor for the time step "dt"
  double& Dynamics::timeStep() {return timeStep_;}
  ///
  ///Read-only accessor for the time step "dt"
  double const& Dynamics::timeStep() const {return timeStep_;}
  ///
  ///Read-write accessor for the number of steps of the simulation
  int& Dynamics::nbOfSteps() {return nbOfSteps_;}
	///
  ///Read-only accessor for the number of steps of the simulation
  int const& Dynamics::nbOfSteps() const {return nbOfSteps_;}
  ///
  ///Returns the total time "t_f" of the simulation
  double Dynamics::finalTime() const {return timeStep_ * nbOfSteps_;}
  ///
  ///Read-write accessor for the number of steps of the thermalization
  int& Dynamics::thermalizationNbOfSteps() {return thermalizationNbOfSteps_;}
  ///
  ///Read-only accessor for the number of steps of the thermalization
  int const& Dynamics::thermalizationNbOfSteps() const {return thermalizationNbOfSteps_;}
  ///
  ///Read-write accessor for the number of steps of the burnIn
  int& Dynamics::burninNbOfSteps() {return burninNbOfSteps_;}
  ///
  ///Read-only accessor for the number of steps of the burnIn
  int const& Dynamics::burninNbOfSteps() const {return burninNbOfSteps_;}

  const std::shared_ptr<RNG>& Dynamics::rng() const {return rng_;}

  std::shared_ptr<RNG>& Dynamics::rng() {return rng_;}
  
  ///
  ///Returns the temperature
  double const& Dynamics::temperature() const {return temperature_;}
  ///
  ///Read-only access for the temperature
  double const& Dynamics::temperatureLeft() const {return temperature_;}
  ///
  ///Read-only access for the temperature
  double const& Dynamics::temperatureRight() const {return temperature_;}
  ///Returns the difference between the mean temperature and the one at the left end
  ///This is equal to {eta}
  double Dynamics::deltaTemperature() const
  {
    return (temperatureLeft() - temperatureRight())/2;
  }
  ///
  ///Read-only access for the inverse temperature
  double const& Dynamics::beta() const {return beta_;}
  ///
  ///Read-only access for the inverse temperature
  double const& Dynamics::betaLeft() const {return beta_;}
  ///
  ///Read-only access for the inverse temperature
  double const& Dynamics::betaRight() const {return beta_;}

  /*///
  ///Read-only accessor for the external force
  Vector<double> const& Dynamics::externalForce() const {return externalForce_;}
  ///
  ///Write-read accessor for the external force
  Vector<double>& Dynamics::externalForce(){return externalForce_;}*/
  /*///
  ///Read-only accessor for the i-th component of the external force
  double const& Dynamics::externalForce(int const& i) const {return externalForce_(i);}
  ///
  ///Write-read accessor for the i-th component of the external force
  double& Dynamics::externalForce(int const& i) {return externalForce_(i);}*/
	///
	///Returns the pointer of the Galerkin instance of *this
  Galerkin* Dynamics::galerkin() {return galerkin_;}




  /*///
  ///Set the force and the potential energy of a particle to zero
  void Dynamics::resetForce(Particle& particle) const
  {
    particle.potentialEnergy() = 0;
    particle.force() = externalForce();
    particle.virial() = 0.;
  }*/


	///
  ///Standard first part of the numerical integration : half upadate on "p" and update on "q"
  void Dynamics::verletFirstPart(Particle& particle)
  {
    particle.momentum() += timeStep_ * particle.force() / 2;
    particle.position() += timeStep_ * particle.momentum() / particle.mass();
  }

  ///
  ///Standard second part of the numerical integration : half upadate on "p"
  void Dynamics::verletSecondPart(Particle& particle)
  {
    particle.momentum() += timeStep_ * particle.force() / 2;
  }

  ///
  ///Standard first part of the numerical integration : half upadate on "p" and update on "q"
  void Dynamics::updateBefore(Particle& particle)
  {
    particle.momentum() += timeStep_ * particle.force() / 2;
    particle.position() += timeStep_ * particle.momentum() / particle.mass();
  }

  ///
  ///Standard second part of the numerical integration : half upadate on "p"
  void Dynamics::updateAfter(Particle& particle)
  {
    particle.momentum() += timeStep_ * particle.force() / 2;
  }



}
#endif
