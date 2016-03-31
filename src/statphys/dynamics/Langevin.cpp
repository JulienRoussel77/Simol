#include "Langevin.hpp"

namespace simol
{
  ///Constructor for the Langevin Dynamics with thermostats everywhere
  Langevin::Langevin(Input const& input):
    UniformStochasticDynamics(input),
    gamma_(input.gamma())
  {
    galerkin_ = createLangevinGalerkin(input);
  }

  ///
  ///Print de name of the dynamics
  void Langevin::printName() const
  {
    std::cout << "Dynamics = Langevin" << std::endl;
  }

  ///
  ///Read-only accessor for the intensity of the Orstein-Uhlenbeck process
  const double& Langevin::gamma() const {return gamma_;}

  ///
  ///Returns the amplitude of the brownian motion
  double Langevin::sigma() const
  {
    return sqrt(2*gamma_ / beta_);
  }

	///After refers to the fact that this step comes after updating the forces
	///Proceeds to the second half of the Verlet scheme, then integrate the O-U analytically
  void Langevin::updateAfter(Particle& particle)
  {
    particle.momentum() += timeStep_ * particle.force() / 2;
    double alpha = exp(- gamma_ / particle.mass() * timeStep_);
    particle.momentum() = alpha * particle.momentum() + sqrt((1-pow(alpha, 2))/beta_*particle.mass()) * rng_->gaussian();
  }

}
