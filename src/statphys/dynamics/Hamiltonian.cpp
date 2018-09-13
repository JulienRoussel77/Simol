
#include "simol/statphys/dynamics/Hamiltonian.hpp"

namespace simol
{
  //! Constructs a Hamiltonian dynamics (constant energy)
  Hamiltonian::Hamiltonian(Input const& input):
    Dynamics(input)
  {}
  
  void Hamiltonian::simulate(System& syst) const
  {    
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      verletFirstPart(syst(iOfParticle));
    syst.computeAllForces();

    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      verletSecondPart(syst(iOfParticle));
  }

}
