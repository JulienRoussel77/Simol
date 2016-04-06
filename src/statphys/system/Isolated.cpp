#include "Isolated.hpp"

namespace simol
{
  
  //### Isolated ###
  
  Isolated::Isolated(Input const& input):
    System(input)
  {
    getParticle() = Particle(input.mass(), input.initialPosition(), input.initialMomentum());
    cout << getParticle().mass() << endl;
    cout << getParticle().mass() << "  " << getParticle().position() << "  " << getParticle().momentum() << endl;
  }
  
  void Isolated::printName() const
  {
    cout << "System = Isolated" << endl;
  }
      
  
  void Isolated::computeAllForces(Dynamics const& dyna)
  {
    dyna.resetForce(getParticle());
    getParticle().potentialEnergy() = potential(getParticle().position());
    getParticle().force() = force(getParticle().position());
  }
  
  //void Isolated::computeFinalOutput(Output& /*output*/, Dynamics const& /*dyna*/)
  //{}
  
}