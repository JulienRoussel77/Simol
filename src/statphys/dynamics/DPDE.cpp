#include "DPDE.hpp"

#include "UniformStochasticDynamics.hpp"

using std::cout;
using std::endl;

namespace simol
{
  DPDE::DPDE(Input const&  input):
    UniformStochasticDynamics(input),
    heatCapacity_(input.heatCapacity()),
    gamma_(input.gamma())
  {}

  double& DPDE::gamma()
  {
    return gamma_;
  }

  double& DPDE::heatCapacity()
  {
    return heatCapacity_;
  }

  double DPDE::gamma_DPDE(double intEnergy)
  {
    return gamma_ * heatCapacity()*temperature()/intEnergy;
  }

  double DPDE::sigma() const
  {
    return sqrt(2*gamma_*temperature());
    // TO DO : ce serait mieux d'avoir return alpha_ pour ne pas recalculer cela a tous les pas de temps !
  }

  void DPDE::energyReinjection(Particle& particle)
  {
    double old_kin_energy = particle.kineticEnergy();
    double local_gamma_DPDE = gamma_DPDE(particle.internalEnergy());
    double alpha = exp(- local_gamma_DPDE / particle.mass() * timeStep_);
    particle.momentum() = alpha * particle.momentum() + sqrt((1-pow(alpha, 2))*particle.mass()*temperature()*gamma()/local_gamma_DPDE) * rng_->gaussian();
    // TO DO : ne pas recalculer la magnitude du terme de fluctuation a chaque pas...
    double new_kin_energy = particle.kineticEnergy();
    particle.internalEnergy() += old_kin_energy-new_kin_energy;
  }

}
