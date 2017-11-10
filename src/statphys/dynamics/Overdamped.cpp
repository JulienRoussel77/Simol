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
  ///Proceeds to an Euler-Maruyama scheme (in comment the second order is available)
  void Overdamped::updatePosition(Particle& particle) const
  {
    //particle.position() += timeStep_ * particle.force() + sqrt(2 * timeStep_ / beta_) * rng_->gaussian();
    DVec newGaussian = rng_->gaussian();
    particle.position() += timeStep_ * particle.force() + sqrt(timeStep_ / (2*parameters_.beta())) * (newGaussian + particle.oldGaussian());
    particle.oldGaussian() = newGaussian;
  }
  
  void Overdamped::simulate(System& syst) const
  {
    //cout << "Overdamped::simulate( System& syst)" << endl;
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
    {
      //cout << iOfParticle << " " << syst(iOfParticle).position().adjoint() << " " << syst(iOfParticle).momentum().adjoint() << " " << syst(iOfParticle).force().adjoint() << endl;
      updatePosition(syst(iOfParticle));
      //cout << iOfParticle << " " << syst(iOfParticle).position().adjoint() << " " << syst(iOfParticle).momentum().adjoint() << " " << syst(iOfParticle).force().adjoint() << endl;
    }
    syst.computeAllForces();
  }
  
  void Overdamped::getThermo(Output& output) const
  {    
    double kineticEnergy = output.dimension() * output.nbOfParticles() / (2 * parameters_.beta());
    
    output.temperature() = 1 / parameters_.beta();
    output.obsTotalEnergy().currentValue() = kineticEnergy + output.potentialEnergy();
  }
  
  ///
  ///Computes the pressure from the kineticEnergy and the totalVirial, these fields must be updated
  void Overdamped::getPressure(Output& output) const
  {
    double kineticEnergy = output.dimension() * output.nbOfParticles() / (2 * parameters_.beta());
    output.pressure() = (2 * kineticEnergy + output.totalVirial()) / (output.dimension() * output.nbOfParticles() * pow(output.latticeParameter(), output.dimension()));
  }

}
