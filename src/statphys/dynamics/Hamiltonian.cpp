
#include "simol/statphys/dynamics/Hamiltonian.hpp"

namespace simol
{
  //! Constructs a Hamiltonian dynamics (constant energy)
  Hamiltonian::Hamiltonian(Input const& input):
    Dynamics(input)
  {}
  
  /*void Hamiltonian::computeThermo(Output& output) const
  {    
    output.temperature() = 2 * output.kineticEnergy() / (output.dimension() * output.nbOfParticles());
    output.totalEnergy() = output.kineticEnergy() + output.potentialEnergy();
    output.pressure() = (2 * output.kineticEnergy() + output.totalVirial()) / (output.dimension() * output.nbOfParticles() * pow(output.latticeParameter(), output.dimension()));
  }*/
  
  void Hamiltonian::simulate(System& syst) const
  {    
    //cout << "Hamiltonian::simulate(System& syst)" << endl;
    //for (ParticleIterator it = syst.begin(); it != syst.end(); ++it)
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      verletFirstPart(syst(iOfParticle));
    syst.computeAllForces();

    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      verletSecondPart(syst(iOfParticle));
  }

}
