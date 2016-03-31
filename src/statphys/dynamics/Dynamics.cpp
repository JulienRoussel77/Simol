#ifndef SIMOL_DYNAMICS_IPP
#define SIMOL_DYNAMICS_IPP

#include "Dynamics.hpp"

using std::cout;
using std::endl;

namespace simol
{



  ///
  ///Constructor
  Dynamics::Dynamics(Input const& input):
    timeStep_(input.timeStep()),
    nbOfIterations_(input.nbOfIterations()),
    nbOfThermalIterations_(input.nbOfThermalIterations()),
    nbOfBurningIterations_(input.nbOfBurningIterations()),
    externalForce_(input.dimension(), 0)
  {
    externalForce_(0) = input.externalForce();

    cout << "externalForce = " << externalForce_(0) << endl;
    cout << "timeStep = " << timeStep() << endl;
    if (nbOfIterations_ < 1e6)
      cout << "nbOfIterations_ = " << nbOfIterations_ << endl;
    else
      cout << "nbOfIterations_ = " << nbOfIterations_/1e6 << "e6" << endl;
    cout << "finalTime = " << nbOfIterations_ * timeStep() << endl;
  }

  void Dynamics::printName() const
  {
    std::cout << "Dynamics = Dynamics" << std::endl;
  }




  ///
  ///Read-write accessor for the time step "dt"
  double& Dynamics::timeStep() {return timeStep_;}
  ///
  ///Read-only accessor for the time step "dt"
  const double& Dynamics::timeStep() const {return timeStep_;}
  ///
  ///Read-write accessor for the number of iterations of the simulation
  size_t& Dynamics::nbOfIterations() {return nbOfIterations_;}
	///
  ///Read-only accessor for the number of iterations of the simulation
  const size_t& Dynamics::nbOfIterations() const {return nbOfIterations_;}
  ///
  ///Returns the total time "t_f" of the simulation
  double Dynamics::finalTime() const {return timeStep_ * nbOfIterations_;}
  ///
  ///Read-write accessor for the number of iterations of the thermalization
  size_t& Dynamics::nbOfThermalIterations() {return nbOfThermalIterations_;}
  ///
  ///Read-only accessor for the number of iterations of the thermalization
  const size_t& Dynamics::nbOfThermalIterations() const {return nbOfThermalIterations_;}
  ///
  ///Read-write accessor for the number of iterations of the burning
  size_t& Dynamics::nbOfBurningIterations() {return nbOfBurningIterations_;}
  ///
  ///Read-only accessor for the number of iterations of the burning
  const size_t& Dynamics::nbOfBurningIterations() const {return nbOfBurningIterations_;}

  const std::shared_ptr<RNG>& Dynamics::rng() const {return rng_;}

  std::shared_ptr<RNG>& Dynamics::rng() {return rng_;}

  ///
  ///Read-only accessor for the external force
  const Vector<double>& Dynamics::externalForce() const {return externalForce_;}
  ///
  ///Write-read accessor for the external force
  Vector<double>& Dynamics::externalForce(){return externalForce_;}
  ///
  ///Read-only accessor for the i-th component of the external force
  const double& Dynamics::externalForce(int const& i) const {return externalForce_(i);}
  ///
  ///Write-read accessor for the i-th component of the external force
  double& Dynamics::externalForce(int const& i) {return externalForce_(i);}
	///
	///Returns the pointer of the Galerkin instance of *this
  Galerkin* Dynamics::galerkin() {return galerkin_;}




  ///
  ///Set the force and the potential energy of a particle to zero
  void Dynamics::resetForce(Particle& particle) const
  {
    particle.potentialEnergy() = 0;
    particle.force() = externalForce();
  }


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

  ///
  ///Analytical integration of an Orstein-Uhlenbeck process of inverse T "localBeta"
  void Dynamics::updateOrsteinUhlenbeck(Particle& particle, double localBeta)
  {
    double alpha = exp(- gamma() / particle.mass() * timeStep_);
    particle.momentum() = alpha * particle.momentum() + sqrt((1-pow(alpha, 2))/localBeta*particle.mass()) * rng_->gaussian();
  }

}
#endif
