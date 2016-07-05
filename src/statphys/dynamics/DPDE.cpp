#include "simol/statphys/dynamics/DPDE.hpp"

#include "simol/statphys/dynamics/LangevinBase.hpp"

using std::cout;
using std::endl;

namespace simol
{
  DPDE::DPDE(Input const&  input):
    LangevinBase(input),
    heatCapacity_(input.heatCapacity())
  {}

  void DPDE::printName() const
  {
    std::cout << "DynamicsType = DPDE" << std::endl;
  }

  double& DPDE::heatCapacity()
  {
    return heatCapacity_;
  }

  double DPDE::gamma_DPDE(double intEnergy)
  {
    return gamma_ * heatCapacity() * temperature() / intEnergy;
  }

  double DPDE::sigma() const
  {
    return sqrt(2 * gamma_ * temperature());
  }

  void DPDE::energyReinjection(Particle& particle)
  {
    double old_kin_energy = particle.kineticEnergy();
    double local_gamma_DPDE = gamma_DPDE(particle.internalEnergy());
    double alpha = exp(- local_gamma_DPDE / particle.mass() * timeStep_);
    particle.momentum() = alpha * particle.momentum() + sqrt((1 - pow(alpha, 2)) * particle.mass() * temperature() * gamma() / local_gamma_DPDE) * rng_->gaussian();
    double new_kin_energy = particle.kineticEnergy();
    particle.internalEnergy() += old_kin_energy - new_kin_energy;
  }

  void DPDE::metropolizedEnergyReinjection(Particle& particle)
  {
    //-- keep previous configuration --
    Vector<double> old_momentum = particle.momentum();
    double old_kin_energy = particle.kineticEnergy();
    //-- propose a new move --
    particle.momentum() += rng_->gaussian();
    //-- compute the Metropolis rate --
    double rate = particle.kineticEnergy() - old_kin_energy;
    //-- acceptance/rejection procedure --
    if (rng_->scalarUniform() > rate)
      particle.momentum() = old_momentum;
  }

}
