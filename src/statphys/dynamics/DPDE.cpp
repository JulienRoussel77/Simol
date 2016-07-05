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

  const double& DPDE::rejectionCount() const
  {
    return rejectionCount_;
  }
  
  double& DPDE::rejectionCount()
  {
    return rejectionCount_;
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
    Vector<double> drift = -pow(sigma(),2)* heatCapacity()/(2*particle.mass())*particle.momentum()/(e0-particle.kineticEnergy());
    return drift;
  }

  void DPDE::incrementRejection()
  {
    rejectionCount_ += 1;
  }

  void DPDE::metropolizedEnergyReinjection(Particle& particle)
  {
    //-- keep previous configuration --
    Vector<double> old_momentum = particle.momentum();
    double old_kin_energy = particle.kineticEnergy();
    double E0 = old_kin_energy + particle.internalEnergy(); 
    //-- propose a new move --
    Vector<double> G = rng_->gaussian();
    Vector<double> b = effectiveDrift(particle);
    particle.momentum() += b*timeStep_ + sigma()*sqrt(timeStep_)*G;
    //-- compute the Metropolis rate --
    double new_kin_energy = particle.kineticEnergy();
    double rate = 0; 
    if (new_kin_energy < E0)
      {
     	rate = particle.mass()*heatCapacity()*(log(E0-new_kin_energy)-log(E0-old_kin_energy)) + 0.5*pow(G.norm(), 2); 
	Vector<double> GG = old_momentum - particle.momentum() - effectiveDrift(particle)*timeStep_;
	rate -= 0.5*pow(GG.norm(), 2)/(timeStep_*pow(sigma(),2));
	rate = exp(rate);
      }
    //-- acceptance/rejection procedure --
    double U = rng_->scalarUniform();
    if (U > rate)
      {
	//-- reject the move --
	incrementRejection();
	//cout << rate << "  " << rejectionCount() << endl;
	//cout << rate << " ; " << E0 << " " << particle.momentum() << " , drift " << b(0) << ", previously " << old_momentum << ", internal = " << particle.internalEnergy() << endl;
	particle.momentum() = old_momentum;
      }
    else 
      {
	//-- update internal energies --
	particle.internalEnergy() = E0 - new_kin_energy;
      }
  }

}
