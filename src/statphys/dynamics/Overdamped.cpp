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
  
  void Overdamped::getThermo(Output& output) const
  {    
    double kineticEnergy = output.dimension() * output.nbOfParticles() / (2 * beta());
    
    output.temperature() = 1 / beta();
    output.totalEnergy() = kineticEnergy + output.potentialEnergy();
  }
  
  ///
  ///Computes the pressure from the kineticEnergy and the totalVirial, these fields must be updated
  void Overdamped::getPressure(Output& output) const
  {
    double kineticEnergy = output.dimension() * output.nbOfParticles() / (2 * beta());
    output.pressure() = (2 * kineticEnergy + output.totalVirial()) / (output.dimension() * output.nbOfParticles() * pow(output.latticeParameter(), output.dimension()));
  }

  
  ///Applies the generator of this dynamics to the basis functions of the CV
  ///Evaluate at the current state of "conifguration"
  void Overdamped::computeGeneratorOnBasis(CVBasis& cvBasis, System const& syst) const
  {
    //Vector<double> result = Vector<double>::Zero(nbOfFunctions());
    for (int iOfFunction = 0; iOfFunction < cvBasis.nbOfFunctions(); iOfFunction++)
      for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
        cvBasis.generatorOnBasisValues_(iOfFunction) += cvBasis.basis_->laplacianQ(syst, iOfParticle, iOfFunction) / beta()
                               + dot(syst(iOfParticle).force(), cvBasis.basis_->gradientQ(syst, iOfParticle, iOfFunction));
  }

}
