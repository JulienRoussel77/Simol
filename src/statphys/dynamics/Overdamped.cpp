#include "Overdamped.hpp"

namespace simol
{
  //#### Overdamped ####

  ///
  ///Constructor for the Overdamped Langevin Dynamics
  Overdamped::Overdamped(Input const& input):
    UniformStochasticDynamics(input)
  {}

  ///Before refers to the fact that this step comes before the forces update
  ///Do nothing
  void Overdamped::updateBefore(Particle& /*particle*/) {}

  // /!\ a changer dans le cas N != 1
  ///After refers to the fact that this step comes after the forces update
  ///Proceeds to a Euler-Maruyama scheme (in comment the second order is available)
  void Overdamped::updateAfter(Particle& particle)
  {
    particle.position() += timeStep_ * particle.force() + sqrt(2*timeStep_/beta_) * rng_->gaussian();

    //assert(particle.force()(0) == force(particle.position())(0));
    /*Vector<double> randomTerm = sqrt(2*timeStep_/beta_) * rng_->gaussian();
    Vector<double> qtilde = particle.position() + .5 * timeStep_ * particle.force() + .5 * randomTerm;
    particle.position() += timeStep_ * force(qtilde) + randomTerm;*/
  }

}
