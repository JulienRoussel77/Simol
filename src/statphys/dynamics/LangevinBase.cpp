#include "simol/statphys/dynamics/LangevinBase.hpp"

namespace simol
{
  ///Constructs a purely virtual Dynamics for Dynamics using RNGs
  LangevinBase::LangevinBase(Input const& input):
    Dynamics(input),
    gamma_(input.gamma()),
    xi_(input.xi())
  {}

  ///
  ///Read-only accessor of the intensity of the O-U process
  const double& LangevinBase::gamma() const {return gamma_;}

  ///Read-only accessor for xi
  const double& LangevinBase::xi() const {return xi_;}

  ///Read-write accessor for xi
  double& LangevinBase::xi() {return xi_;}
  ///
  ///Returns the mean number of steps between 2 random events
  int LangevinBase::xiNbOfSteps()
  {return 1 / (xi_ * timeStep_);}
  ///
  ///Returns true if the dynamics involves a Poisson process (momenta exchange)
  bool LangevinBase::doMomentaExchange() const
  {return xi_ > 0;}
  ///
  ///If the momenta exchange is activated, the times of future events are drawn
  void LangevinBase::initializeCountdown(Particle& particle)
  {
    if (doMomentaExchange())
      particle.countdown() = rng_->scalarExponential() * xiNbOfSteps(); // / (xi() * timestep());
    else
      particle.countdown() = -1;
  }
  ///
  ///Exchanges the momenta of the 2 particles if the time has come
  void LangevinBase::updateMomentaExchange(Particle& particle1, Particle& particle2)
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
  void LangevinBase::updateOrsteinUhlenbeck(Particle& particle, double localBeta)
  {
    double alpha = exp(- gamma() / particle.mass() * timeStep_);
    particle.momentum() = alpha * particle.momentum() + sqrt((1 - pow(alpha, 2)) / localBeta * particle.mass()) * rng_->gaussian();
  }
  
  void LangevinBase::getThermo(Output& output) const
  {
    output.temperature() = 2 * output.kineticEnergy() / (output.dimension() * output.nbOfParticles());
    output.totalEnergy() = output.kineticEnergy() + output.potentialEnergy();
  }
  
  ///
  ///Computes the pressure from the kineticEnergy and the totalVirial, these fields must be updated
  void LangevinBase::getPressure(Output& output) const
  {
    output.pressure() = (2 * output.kineticEnergy() + output.totalVirial()) / (output.dimension() * output.nbOfParticles() * pow(output.latticeParameter(), output.dimension()));
  }


}
