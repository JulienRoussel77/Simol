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

  Vector<double> DPDE::effectiveDrift(Particle& particle)
  {
    double e0 = particle.totalEnergyDPDE(); 
    Vector<double> drift = -pow(sigma(),2)/(2*particle.mass())*particle.momentum()/(e0-particle.kineticEnergy());
    return drift;
  }

  void DPDE::metropolizedEnergyReinjection(Particle& particle)
  {
    //-- keep previous configuration --
    Vector<double> old_momentum = particle.momentum();
    double E0 = particle.totalEnergyDPDE(); 
    //-- propose a new move --
    Vector<double> G = rng_->gaussian();
    particle.momentum() += effectiveDrift(particle)*timeStep_ + sigma()*sqrt(timeStep_)*G;
    //-- compute the Metropolis rate --
    double new_kin_energy = particle.kineticEnergy();
    double rate = 0; 
    if (new_kin_energy < E0)
      {
	rate = 1; //pow(G.norm(), 2); 
      }
    //-- acceptance/rejection procedure --
    if (rng_->scalarUniform() < rate)
      particle.momentum() = old_momentum;
  }

}
