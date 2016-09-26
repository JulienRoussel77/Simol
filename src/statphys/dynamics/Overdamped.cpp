#include "simol/statphys/dynamics/Overdamped.hpp"

namespace simol
{
  //#### Overdamped ####

  ///
  ///Constructor for the Overdamped Langevin Dynamics
  Overdamped::Overdamped(Input const& input):
    Dynamics(input)
  {}

  // /!\ a changer dans le cas N != 1
  ///After refers to the fact that this step comes after the forces update
  ///Proceeds to a Euler-Maruyama scheme (in comment the second order is available)
  void Overdamped::updatePosition(Particle& particle)
  {
    particle.position() += timeStep_ * particle.force() + sqrt(2 * timeStep_ / beta_) * rng_->gaussian();
  }

}
