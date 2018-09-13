#include "simol/dynamics/LangevinBase.hpp"

namespace simol
{
  ///Constructs a purely virtual Dynamics for Dynamics using RNGs
  LangevinBase::LangevinBase(Input const& input):
    Dynamics(input)
  {}

  ///
  ///Read-only accessor of the intensity of the O-U process
  const double& LangevinBase::gamma() const {return parameters_.gamma();}

  ///Read-only accessor for xi
  const double& LangevinBase::xi() const {return parameters_.xi();}
  ///
  ///Returns the mean number of steps between 2 random events
  int LangevinBase::xiNbOfSteps() const
  {return 1 / (xi() * timeStep_);}
  ///
  ///Returns true if the dynamics involves a Poisson process (momenta exchange)
  bool LangevinBase::doMomentaExchange() const
  {return xi() > 0;}
  ///
  ///If the momenta exchange is activated, the times of future events are drawn
  void LangevinBase::initializeCountdown(Particle& particle) const
  {
    if (doMomentaExchange())
      particle.countdown() = rng_->scalarExponential() * xiNbOfSteps(); // / (xi() * timestep());
    else
      particle.countdown() = -1;
  }
  ///
  ///Exchanges the momenta of the 2 particles if the time has come
  void LangevinBase::updateMomentaExchange(Particle& particle1, Particle& particle2) const
  {
    if (particle2.countdown() == 0)
    {
      DVec temp = particle2.momentum();
      particle2.momentum() = particle1.momentum();
      particle1.momentum() = temp;
      particle2.countdown() = rng_->scalarExponential() * xiNbOfSteps(); // / (xi() * timestep());
    }
    else
      particle2.countdown()--;
  }

  ///
  ///Analytical integration of an Orstein-Uhlenbeck process of inverse T "localBeta"
  void LangevinBase::updateOrsteinUhlenbeck(Particle& particle, double localBeta, double localTimeStep) const
  {
    double alpha = exp(- parameters_.gamma() / particle.mass() * localTimeStep);
    particle.momentum() = alpha * particle.momentum() + sqrt((1 - pow(alpha, 2)) / localBeta * particle.mass()) * rng_->gaussian();
  }
  
  void LangevinBase::simulate(System& syst) const
  {
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      verletFirstPart(syst(iOfParticle));
    syst.computeAllForces();
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      verletSecondPart(syst(iOfParticle));
  }
  



}
