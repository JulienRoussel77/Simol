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
  
  ///Applies the generator of this dynamics to the basis functions of the CV
  ///Evaluate at the current state of "conifguration"
  void Overdamped::computeGeneratorOnBasis(CVBasis& cvBasis, vector<Particle> const& configuration) const
  {
    int nbOfParticles = (int)configuration.size();
    //Vector<double> result = Vector<double>::Zero(nbOfFunctions());
    for (int iOfFunction = 0; iOfFunction < cvBasis.nbOfFunctions(); iOfFunction++)
      for (int iOfParticle = 0; iOfParticle < nbOfParticles; iOfParticle++)
        cvBasis.generatorOnBasisValues_(iOfFunction) += cvBasis.basis_->laplacianQ(configuration, iOfParticle, iOfFunction) / beta()
                               + dot(configuration[iOfParticle].force(), cvBasis.basis_->gradientQ(configuration, iOfParticle, iOfFunction));
  }

}
