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
  
  void Langevin::simulate(System& syst) const
  {
    
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      verletFirstPart(syst(iOfParticle));
    syst.computeAllForces();

    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      verletSecondPart(syst(iOfParticle));

    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      updateOrsteinUhlenbeck(syst(iOfParticle), beta(), timeStep());
  }
  
  ///Constructor for the Constrained Langevin Dynamics with thermostats everywhere
  ConstrainedLangevin::ConstrainedLangevin(Input const& input):
    Langevin(input)
  {}

void ConstrainedLangevin::simulate(System& syst) const
{
  double alpha = exp(- gamma() / syst(0).mass() * timeStep()/2);
  syst.lagrangeMultiplier() = 2 * (1-alpha) * drift();
  
  for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
    updateOrsteinUhlenbeck(syst(iOfParticle), beta(), timeStep()/2);
  
  // Becareful here we assume that all the particles share the same mass !
  // We analiticaly determine the non-martingale part of the Lagrange multiplier for the O-U process
  
  //syst.lagrangeMultiplier() += (1-alpha) * drift();
  
  syst.enforceConstraint(drift(), parameters(), false);
  //syst.enforceConstraint(syst.lagrangeMultiplier(), drift());
  
  for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
    updateMomentum(syst(iOfParticle));
  syst.enforceConstraint(drift(), parameters(), true);
  
  for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
    updatePosition(syst(iOfParticle));
  syst.computeAllForces();

  for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
    verletSecondPart(syst(iOfParticle));
  syst.enforceConstraint(drift(), parameters(), true);
  
  for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
    updateOrsteinUhlenbeck(syst(iOfParticle), beta(), timeStep()/2);
  
  //syst.lagrangeMultiplier() += (1-alpha) * drift();
  syst.enforceConstraint(drift(), parameters(), false);
  //syst.enforceConstraint(syst.lagrangeMultiplier(), drift());
  
  syst.lagrangeMultiplier() /= timeStep();
}

}
