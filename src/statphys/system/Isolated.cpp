#include "simol/statphys/system/Isolated.hpp"

namespace simol
{

  Isolated::Isolated(Input const& input):
    System(input)
  {
    configuration_[0] = new Particle(input.mass(), input.initialPosition(), input.initialMomentum());
  }

  void Isolated::printName() const
  {
    cout << "System = Isolated" << endl;
  }
  
  void Isolated::samplePositions(DynamicsParameters const& dynaPara)
  {
    cout << " - Sampling the positions..." << endl;
    getParticle(0).position(0) = drawPotLaw(dynaPara.beta());
  }

  void Isolated::computeAllForces()
  {
    getParticle().resetForce(potential());
    getParticle().potentialEnergy() = potential(getParticle().position());
    getParticle().force() = totalForce(getParticle().position());
  }
  
  /*void Isolated::computeKineticEnergy(Output& output) const
  {
    output.kineticEnergy() = getParticle().kineticEnergy();
  }
  
  void Isolated::computePotentialEnergy(Output& output) const
  {
    output.potentialEnergy() = getParticle().potentialEnergy();
  }*/
  
  /*void Isolated::computePressure(Output& output) const
  {
    output.totalVirial() = getParticle().virial();
  }*/
  
  


}
