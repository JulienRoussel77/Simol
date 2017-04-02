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
    cout << "--> q = " << getParticle(0).position(0) << endl;
  }

  void Isolated::computeAllForces()
  {
    getParticle().resetForce(potential());
    getParticle().potentialEnergy() = potential(getParticle().position());
    getParticle().force() = totalForce(getParticle().position());
  }
  
  
  


}
