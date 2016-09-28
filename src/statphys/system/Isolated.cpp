#include "simol/statphys/system/Isolated.hpp"

namespace simol
{

  Isolated::Isolated(Input const& input):
    System(input)
  {
    getParticle() = Particle(input.mass(), input.initialPosition(), input.initialMomentum());
    //cout << getParticle().mass() << endl;
    //cout << getParticle().mass() << "  " << getParticle().position() << "  " << getParticle().momentum() << endl;
  }

  void Isolated::printName() const
  {
    cout << "System = Isolated" << endl;
  }

  void Isolated::computeAllForces()
  {
    getParticle().resetForce(potential());
    getParticle().potentialEnergy() = potential(getParticle().position());
    getParticle().force() = totalForce(getParticle().position());
  }
  
  void Isolated::computeKineticEnergy(Output& output) const
  {
    output.kineticEnergy() = getParticle().kineticEnergy();
  }
  
  void Isolated::computePotentialEnergy(Output& output) const
  {
    output.potentialEnergy() = getParticle().potentialEnergy();
  }
  
  /*void Isolated::computePressure(Output& output) const
  {
    output.totalVirial() = getParticle().virial();
  }*/
  
  
  void Isolated::thermalize(Dynamics& /*model*/)
  {}

}
