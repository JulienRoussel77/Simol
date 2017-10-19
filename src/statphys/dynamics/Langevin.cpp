#include "simol/statphys/dynamics/Langevin.hpp"

namespace simol
{
  ///Constructor for the Langevin Dynamics with thermostats everywhere
  Langevin::Langevin(Input const& input):
    LangevinBase(input)
  {}


  ///
  ///Returns the amplitude of the brownian motion
  double Langevin::sigma() const
  {
    return sqrt(2 * parameters_.gamma() / parameters_.beta());
  }

  /*///After refers to the fact that this step comes after updating the forces
  ///Proceeds to the second half of the Verlet scheme, then integrate the O-U analytically
  void Langevin::updateAfter(Particle& particle)
  {
    particle.momentum() += timeStep_ * particle.force() / 2;
    double alpha = exp(- gamma_ / particle.mass() * timeStep_);
    particle.momentum() = alpha * particle.momentum() + sqrt((1 - pow(alpha, 2)) / beta_ * particle.mass()) * rng_->gaussian();
  }*/
  
  
  
  ///Constructor for the Constrained Langevin Dynamics with thermostats everywhere
  ConstrainedLangevin::ConstrainedLangevin(Input const& input):
    Langevin(input)
  {}
  
  /*void ConstrainedLangevin::halfUpdateOrsteinUhlenbeck(System& syst) const
  {
    // Becareful here we assume that all the particles share the same mass !
    double alpha = exp(- gamma() / syst(0).mass() * timeStep()/2);
    double localLagrangeMultiplier = 0;
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
    {
      localLagrangeMultiplier += 
      syst(iOfParticle).momentum() = alpha * particle.momentum() + sqrt((1 - pow(alpha, 2)) / localBeta * particle.mass()) * rng_->gaussian();
  }*/
  


}
