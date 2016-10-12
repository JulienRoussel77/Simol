#include "simol/statphys/simulation/Simulation.hpp"

namespace simol
{


  //--- DPDE dynamics ---
  void thermalize(DPDE& dyna, System& syst)
  {
    //-- Verlet part --
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      dyna.verletFirstPart(syst(iOfParticle));
    syst.computeAllForces();
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      dyna.secondPartThermalization(syst(iOfParticle));
  }

  ///
  /// A SUPPRIMER ? NE DEVRAIT PAS DEPENDRE DU SYSTEME...
  void simulate(DPDE& dyna, NBody& syst)
  {
    //-- Verlet part --
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      dyna.verletFirstPart(syst(iOfParticle));
    syst.computeAllForces();
    for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      dyna.verletSecondPart(syst(iOfParticle));
    //-- fluctuation/dissipation --
    dyna.fluctuationDissipation(syst);
    //-- thermal conduction --
    dyna.thermalConduction(syst);
  }
}
